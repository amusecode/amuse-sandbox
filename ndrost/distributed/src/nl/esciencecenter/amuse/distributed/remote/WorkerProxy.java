package nl.esciencecenter.amuse.distributed.remote;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.WorkerDescription;

import ibis.amuse.Daemon;
import ibis.ipl.Ibis;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * In charge of starting a AMUSE "worker" locally, and receiving requests from a WorkerConnection, and forwarding these to the
 * worker via a loopback socket.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerProxy extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(WorkerProxy.class);

    private static final int ACCEPT_TRIES = 20;
    private static final int ACCEPT_TIMEOUT = 1000; // ms

    private static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_" };

    // local socket communication stuff

    private final SocketChannel socket;

    private final Process process;
    private int result = 0;

    private final OutputForwarder out;
    private final OutputForwarder err;

    private final WorkerDescription description;
    private final AmuseConfiguration amuseConfiguration;

    //connection back to AMUSE

    private final Ibis ibis;

    private final File workingDirectory;

    /**
     * List of all hosts used, to give to MPI. Contains duplicates for all machines running multiple worker processes
     */
    private static String[] createHostnameList(WorkerDescription description, IbisIdentifier[] nodes)
            throws DistributedAmuseException {
        String[] hostnames = new String[description.getNrOfWorkers()];
        int next = 0;

        int nrOfProcesses = description.getNrOfWorkers();
        int nrOfNodes = description.getNrOfNodes();

        if (nrOfNodes != nodes.length) {
            throw new DistributedAmuseException("number of nodes used (" + nodes.length
                    + ") not equals to number of nodes required (" + nrOfNodes + ")");
        }

        hostnames = new String[nrOfProcesses];

        if (nodes.length == 1) {
            //just us, job done.
            for (int i = 0; i < nrOfProcesses; i++) {
                hostnames[i] = "localhost";
            }
        } else {
            for (int i = 0; i < nodes.length; i++) {
                String hostname = nodes[i].tagAsString();

                // number of processes per node
                int ppn = nrOfProcesses / nrOfNodes;
                // nrOfWorkers not divideble by nrOfNodes. see if this is
                // a "remainder" node with an extra worker
                if (i < nrOfProcesses % nrOfNodes) {
                    ppn++;
                }
                for (int j = 0; j < ppn; j++) {
                    hostnames[next] = hostname;
                    next++;
                }
            }
            if (next != nrOfProcesses) {
                logger.error("error in setting hostnames. List is of length " + next + " but should be " + nrOfProcesses);
            }
        }
        return hostnames;
    }

    private static Process startWorkerProcess(WorkerDescription description, AmuseConfiguration amuseConfiguration,
            int localSocketPort, String[] hostnames, File workingDirectory) throws Exception {
        File executable;

        if (description.copyWorkerCode()) {
            // executable in workingDirectory
            executable = new File(description.getCodeName());
        } else {
            executable =
                    new File(amuseConfiguration.getAmuseHome() + File.separator + description.getCodeDir() + File.separator
                            + description.getCodeName());
        }

        if (!description.copyWorkerCode() && !executable.canExecute()) {
            throw new DistributedAmuseException(executable + " is not executable, or does not exist");
        }

        File hostFile = File.createTempFile("host", "file").getAbsoluteFile();

        FileWriter hostFileWriter = new FileWriter(hostFile);
        for (String hostname : hostnames) {
            hostFileWriter.write(hostname + "\n");
        }
        hostFileWriter.flush();
        hostFileWriter.close();

        logger.info("host file = " + hostFile);

        ProcessBuilder builder = new ProcessBuilder();

        builder.directory(workingDirectory);

        // make sure there is an "output" directory for a code to put output in
        new File("output").mkdir();

        for (String key : builder.environment().keySet().toArray(new String[0])) {
            for (String blacklistedKey : ENVIRONMENT_BLACKLIST) {
                if (key.startsWith(blacklistedKey)) {
                    builder.environment().remove(key);
                    logger.info("removed " + key + " from environment");
                }
            }
        }
        if (description.getNrOfThreads() > 0) {
            builder.environment().put("OMP_NUM_THREADS", Integer.toString(description.getNrOfThreads()));
        }

        if (!amuseConfiguration.isMpiexecEnabled()) {
            logger.info("not using mpiexec (as it is disabled)");
            if (description.getNrOfWorkers() > 1) {
                throw new DistributedAmuseException("multiple workers (" + description.getNrOfWorkers()
                        + ") requested, but mpiexec disabled in this AMUSE installation");
            }
        } else if (description.getNrOfNodes() == 1) {
            // no need for machine file, set number of processes.
            builder.command(amuseConfiguration.getMpiexec(), "-n", Integer.toString(description.getNrOfWorkers()));
        } else {
            // use machine file
            builder.command(amuseConfiguration.getMpiexec(), "-machinefile", hostFile.getAbsolutePath());
        }

        if (description.copyWorkerCode()) {
            // run executable via amuse.sh script to set python path and interpreter
            builder.command().add(new File(amuseConfiguration.getAmuseHome(), "amuse.sh").getAbsolutePath());
        }

        // executable and port options
        builder.command().add(executable.toString());
        builder.command().add(Integer.toString(localSocketPort));

        logger.info("starting worker process, command = " + builder.command());

        //start process and return
        return builder.start();
    }

    private static SocketChannel acceptConnection(ServerSocketChannel serverSocket) throws IOException {
        serverSocket.configureBlocking(false);
        for (int i = 0; i < ACCEPT_TRIES; i++) {
            SocketChannel result = serverSocket.accept();

            if (result != null) {
                result.socket().setTcpNoDelay(true);
                return result;
            }
            try {
                Thread.sleep(ACCEPT_TIMEOUT);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }
        throw new IOException("worker not started, socket connection failed to initialize");
    }

    /**
     * Starts a worker proxy. Make take a while.
     */
    WorkerProxy(WorkerDescription description, AmuseConfiguration amuseConfiguration, IbisIdentifier[] nodes, Ibis ibis,
            File workingDirectory) throws Exception {
        this.description = description;
        this.amuseConfiguration = amuseConfiguration;
        this.ibis = ibis;
        this.workingDirectory = workingDirectory;

        String[] hostnames = createHostnameList(description, nodes);

        ServerSocketChannel serverSocket = ServerSocketChannel.open();
        serverSocket.bind(new InetSocketAddress(InetAddress.getByName(null), 0));

        //create process
        process =
                startWorkerProcess(description, amuseConfiguration, serverSocket.socket().getLocalPort(), hostnames,
                        workingDirectory);

        //attach streams
        out = new OutputForwarder(process.getInputStream(), description.getStdoutFile(), ibis);
        err = new OutputForwarder(process.getErrorStream(), description.getStderrFile(), ibis);

        logger.info("process started");

        socket = acceptConnection(serverSocket);
        serverSocket.close();

        logger.info("connection with local worker process established");

        //create a connection back to the amuse process via the ibis there.

        //start a thread to start handling amuse requests
        setName("Worker Proxy for " + description.getID());
        setDaemon(true);
        start();
    }

    public synchronized void end() {
        if (process != null) {
            process.destroy();

            try {
                result = process.exitValue();
                logger.info("Process ended with result " + result);
            } catch (IllegalThreadStateException e) {
                logger.error("Process not ended after process.destroy()!");
            }
        }

        if (out != null) {
            // wait for out and err a bit
            try {
                out.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

        if (err != null) {
            try {
                err.join(1000);
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

    }

    /**
     * First connects to the "home" AMUSE ibis. Then continuously receives a message, performs a call, sends a reply.
     */
    public void run() {
        boolean running = true;
        long start, finish;

        AmuseMessage requestMessage = new AmuseMessage();
        AmuseMessage resultMessage = new AmuseMessage();

        SendPort sendPort = null;
        ReceivePort receivePort = null;

        try {
            sendPort = ibis.createSendPort(Daemon.portType);
            receivePort = ibis.createReceivePort(Daemon.portType, description.getID());

            receivePort.enableConnections();

            // get address of amuse node
            logger.debug("waiting for result of amuse election");
            IbisIdentifier amuse = ibis.registry().getElectionResult("amuse");

            logger.debug("connecting to receive port of worker at amuse node");
            sendPort.connect(amuse, description.getID());
            logger.debug("connected, saying hello");
            WriteMessage helloMessage = sendPort.newMessage();
            helloMessage.writeObject(receivePort.identifier());
            helloMessage.finish();

            while (running) {
                start = System.currentTimeMillis();
                logger.debug("Receiving call message");
                ReadMessage readMessage = receivePort.receive();

                if (readMessage == null) {
                    throw new IOException("cannot get request from worker");
                }

                logger.debug("Reading call request from IPL message");

                requestMessage.readFrom(readMessage);

                readMessage.finish();

                int functionID = requestMessage.getFunctionID();

                if (functionID == AmuseMessage.FUNCTION_ID_STOP) {
                    // final request handled
                    running = false;
                }

                // clean message
                resultMessage.clear();

                logger.debug("Performing call for function " + functionID);

                // perform call. Will put result in result message
                logger.debug("performing call with function ID " + requestMessage.getFunctionID());
                requestMessage.writeTo(socket);
                resultMessage.readFrom(socket);
                logger.debug("done performing call with function ID " + requestMessage.getFunctionID() + " error = "
                        + resultMessage.getError());

                logger.debug("result: " + resultMessage);

                WriteMessage writeMessage = sendPort.newMessage();

                resultMessage.writeTo(writeMessage);

                writeMessage.finish();

                logger.debug("Done performing call for function " + functionID);
                finish = System.currentTimeMillis();

                if (logger.isDebugEnabled()) {
                    logger.debug("Call took " + (finish - start) + " ms");
                }

            }
        } catch (Exception e) {
            logger.error("Error while handling request, stopping worker", e);
        }
        if (sendPort != null) {
            try {
                sendPort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }
        if (receivePort != null) {
            try {
                receivePort.close();
            } catch (IOException e) {
                //IGNORE
            }
        }

    }

}
