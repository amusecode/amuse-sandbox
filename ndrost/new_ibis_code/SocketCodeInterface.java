package ibis.amuse;

import ibis.ipl.IbisIdentifier;
import ibis.util.RunProcess;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class representing a code that is used via a loopback socket interface.
 * 
 * @author Niels Drost
 * 
 */
public class SocketCodeInterface implements CodeInterface {

    private static final Logger logger = LoggerFactory.getLogger(SocketCodeInterface.class);

    private static final int ACCEPT_TRIES = 20;
    private static final int ACCEPT_TIMEOUT = 1000; // ms

    private static final String[] ENVIRONMENT_BLACKLIST = { "JOB_ID", "PE_", "PRUN_", "JOB_NAME", "JOB_SCRIPT", "OMPI_" };

    private File executable;

    // local socket communication stuff

    private ServerSocketChannel serverSocket;

    private SocketChannel socket;

    private Process process;
    private int result = 0;

    private AmuseMessage requestMessage;
    private AmuseMessage resultMessage;

    private OutputPrefixForwarder out;
    private OutputPrefixForwarder err;

    private final WorkerInfo info;

    private final String[] hostnames;
    private final int[] ports;

    SocketCodeInterface(WorkerInfo info, IbisIdentifier[] ibises) throws CodeException {
        this.info = info;
        int next = 0;
        int nrOfProcesses = info.getNrOfProcesses();

        hostnames = new String[nrOfProcesses];
        ports = new int[info.getNrOfProcesses()];

        if (ibises == null) {
            for (int i = 0; i < nrOfProcesses; i++) {
                hostnames[i] = "localhost";
                ports[i] = 0;
            }
        } else {
            for (int i = 0; i < ibises.length; i++) {
                String[] hostAndPort = ibises[i].tagAsString().split(",");
                if (hostAndPort.length != 2) {
                    throw new CodeException("could not get host and port from string " + ibises[i].tagAsString());
                }
                String hostname = hostAndPort[0];
                int port = Integer.parseInt(hostAndPort[1]);

                // number of processes per node
                int workerCount = info.getNrOfProcesses() / ibises.length;
                // nrOfWorkers not divideble by number of hosts. see if this is
                // a
                // "remainder" node with an extra worker
                if (i < info.getNrOfProcesses() % ibises.length) {
                    workerCount++;
                }
                for (int j = 0; j < workerCount; j++) {
                    hostnames[next] = hostname;
                    ports[next] = port;

                    next++;
                }
            }
            if (next != info.getNrOfProcesses()) {
                logger.error("error in setting ibises. List is of length " + next + " but should be " + nrOfProcesses);
            }
        }
    }

    public void init() throws CodeException {
        if (info.getMpdboot() != null) {
            RunProcess p = new RunProcess(info.getMpdboot(), "--chkup");
            p.run();
            logger.info("run mpdboot --chkup, result = " + p.getExitStatus());
            logger.info("stdout = " + new String(p.getStdout()));
            logger.info("stderr = " + new String(p.getStderr()));
        }

        executable = new File(info.getAmuseHome() + File.separator + info.getCodeDir() + File.separator
                + info.getCodeName());
        if (!executable.isFile()) {
            throw new CodeException("Cannot find executable for code " + info.getCodeName() + ": " + executable);
        }

        if (!executable.canExecute()) {
            throw new CodeException(executable + " is not executable");
        }

        try {
            serverSocket = ServerSocketChannel.open();
            serverSocket.socket().bind(null);

            File hostFile = File.createTempFile("host", "file").getAbsoluteFile();

            FileWriter hostFileWriter = new FileWriter(hostFile);
            for (String hostname : hostnames) {
                hostFileWriter.write(hostname + "\n");
            }
            hostFileWriter.flush();
            hostFileWriter.close();

            logger.info("host file = " + hostFile);

            File portFile = File.createTempFile("port", "file").getAbsoluteFile();

            FileWriter portFileWriter = new FileWriter(portFile);
            for (int port : ports) {
                portFileWriter.write(port + "\n");
            }
            portFileWriter.flush();
            portFileWriter.close();

            logger.info("port file = " + portFile);

            ProcessBuilder builder = new ProcessBuilder();

            for (String key : builder.environment().keySet().toArray(new String[0])) {
                for (String blacklistedKey : ENVIRONMENT_BLACKLIST) {
                    if (key.startsWith(blacklistedKey)) {
                        builder.environment().remove(key);
                        logger.info("removed " + key + " from environment");
                    }
                }
            }

            builder.environment().put("OMPI_IBIS_PROFILING_PORT_FILE", portFile.getAbsolutePath());

            if (info.getMpiexec() == null || info.getMpiexec().equalsIgnoreCase("none")) {
                logger.info("not using mpiexec (as it is not set, or set to 'none')");
                if (info.getNrOfProcesses() > 1) {
                    logger.warn("this probably won't work, as more than one worker is requested");
                }
                builder.command(executable.toString(), Integer.toString(serverSocket.socket().getLocalPort()));
            } else if (info.getNrOfNodes() == 1) {
                // no need for machine file
                builder.command(info.getMpiexec(), "-n", Integer.toString(info.getNrOfProcesses()),
                        executable.toString(), Integer.toString(serverSocket.socket().getLocalPort()));

            } else {
                builder.command(info.getMpiexec(), "-machinefile", hostFile.getAbsolutePath(), executable.toString(),
                        Integer.toString(serverSocket.socket().getLocalPort()));
            }

            // make sure there is an "output" directory for a code to put output
            // in
            new File("output").mkdir();

            logger.info("starting worker process, command = " + builder.command());

            synchronized (this) {
                process = builder.start();

                out = new OutputPrefixForwarder(process.getInputStream(), System.out, info.getID());
                err = new OutputPrefixForwarder(process.getErrorStream(), System.err, info.getID());

            }

            logger.info("process started");

            socket = acceptConnection(serverSocket);

            logger.info("connection established");

        } catch (IOException e) {
            throw new CodeException("Cannot initialize socket code interface for " + info.getID() + ": " + e.getMessage(), e);
        }
    }

    private static SocketChannel acceptConnection(ServerSocketChannel serverSocket) throws IOException {
        serverSocket.configureBlocking(false);
        for (int i = 0; i < ACCEPT_TRIES; i++) {
            SocketChannel result = serverSocket.accept();

            if (result != null) {
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

    @Override
    public void call() throws CodeException {
        try {
            logger.debug("performing call with function ID " + requestMessage.getFunctionID());
            requestMessage.writeTo(socket);
            resultMessage.readFrom(socket);
            logger.debug("done performing call with function ID " + requestMessage.getFunctionID() + " error = "
                    + resultMessage.getError());
        } catch (IOException e) {
            throw new CodeException("communication with code failed", e);
        }
    }

    @Override
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

    @Override
    public void setRequestMessage(AmuseMessage message) {
        this.requestMessage = message;
    }

    @Override
    public void setResultMessage(AmuseMessage message) {
        this.resultMessage = message;
    }

    @Override
    public synchronized int getResult() {
        return result;
    }

}
