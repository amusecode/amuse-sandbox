/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed.local;

import nl.esciencecenter.amuse.distributed.AmuseMessage;
import nl.esciencecenter.amuse.distributed.WorkerDescription;

import ibis.amuse.Daemon;
import ibis.amuse.Worker;
import ibis.ipl.ConnectionClosedException;
import ibis.ipl.Ibis;
import ibis.ipl.ReadMessage;
import ibis.ipl.ReceivePort;
import ibis.ipl.ReceivePortIdentifier;
import ibis.ipl.ReceiveTimedOutException;
import ibis.ipl.SendPort;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class responsible for taking method invocation messages from AMUSE, and forwarding them to a remote worker proxy.
 * @author Niels Drost
 * 
 */
public class WorkerConnection extends Thread {

    private static final Logger logger = LoggerFactory.getLogger(Worker.class);

    public static final int CONNECT_TIMEOUT = 60000;

    private static int nextID = 0;
    
    private static int getNextID() {
        return nextID++;
    }
    
    private final SocketChannel socket;

    private final String id;

    private final Ibis ibis;

    private final ReceivePort receivePort;

    private final SendPort sendPort;

    private final AmuseMessage initRequest;

    private final AmuseWorkerJob job;

    /*
     * Initializes worker by reading settings from amuse, deploying the worker
     * process on a (possibly remote) machine, and waiting for a connection from
     * the worker
     */
    WorkerConnection(SocketChannel socket, DistributedAmuse distributedAmuse) throws Exception {
        this.socket = socket;
        this.ibis = distributedAmuse.getNetwork().getIbis();

        if (logger.isDebugEnabled()) {
            logger.debug("New worker connection from " + socket.socket().getRemoteSocketAddress());
        }

        // read initialization call

        initRequest = new AmuseMessage();
        initRequest.readFrom(socket);

        if (initRequest.getFunctionID() != AmuseMessage.FUNCTION_ID_INIT) {
            throw new IOException("first call to worker must be init function");
        }

        String codeName = initRequest.getString(0);
        String codeDir = initRequest.getString(1);
        //String hostname = initRequest.getString(2);
        String stdoutFile = initRequest.getString(3);
        String stderrFile = initRequest.getString(4);
        String nodeLabel = initRequest.getString(5);

        int nrOfWorkers = initRequest.getInteger(0);
        int nrOfNodes = initRequest.getInteger(1);
        int nrOfThreads = initRequest.getInteger(2);

        boolean copyWorkerCode = initRequest.getBoolean(0);
        
        // get rid of "ugly" parts of id
        String idName = codeName.replace("_worker", "").replace("_sockets", "");
        
        id = idName + "-" + getNextID();

        //description of the worker, used for both the scheduler and the code proxy to start the worker properly
        WorkerDescription workerDescription = new WorkerDescription(id, codeName, codeDir, stdoutFile, stderrFile, nodeLabel, nrOfWorkers, nrOfNodes, nrOfThreads, copyWorkerCode);
        
        // initialize ibis ports
        receivePort = ibis.createReceivePort(Daemon.portType, id);
        receivePort.enableConnections();
        
        sendPort = ibis.createSendPort(Daemon.portType);

        // start deployment of worker (possibly on remote machine)
        job = distributedAmuse.getScheduler().submitWorkerJob(workerDescription);

        logger.info("New worker submitted: " + this);

        setDaemon(true);
        start();
    }

    public void run() {
        AmuseMessage request = new AmuseMessage();
        AmuseMessage result = new AmuseMessage();
        long start, finish;

        // finish initializing worker
        try {
            
            // wait until job is running
            job.waitUntilStarted();
            
            //read initial "hello" message with identifier
            ReadMessage helloMessage = receivePort.receive(CONNECT_TIMEOUT);
            
            ReceivePortIdentifier remotePort = (ReceivePortIdentifier) helloMessage.readObject();
            
            helloMessage.finish();

            sendPort.connect(remotePort, CONNECT_TIMEOUT, true);

            // do init function at remote worker so it can initialize the code

            // write init message
            WriteMessage initWriteMessage = sendPort.newMessage();
            initRequest.writeTo(initWriteMessage);
            initWriteMessage.finish();

            // read reply
            AmuseMessage initReply = new AmuseMessage();
            ReadMessage initReadMessage = receivePort.receive();
            initReply.readFrom(initReadMessage);
            initReadMessage.finish();

            if (initReply.getError() != null) {
                throw new IOException(initReply.getError());
            }

            // send reply to amuse
            initReply.writeTo(socket);

        } catch (Exception e) {
            if (socket.isOpen() && socket.isConnected()) {
                logger.error("Error on handling call", e);

                // report error to amuse
                AmuseMessage errormessage =
                        new AmuseMessage(initRequest.getCallID(), initRequest.getFunctionID(), initRequest.getCallCount(),
                                "Amuse error: " + e.getMessage());
                try {
                    errormessage.writeTo(socket);
                } catch (IOException e1) {
                    logger.error("Error while returning error message to amuse", e1);
                }
            } else {
                logger.error("Error on handling call, lost connection to AMUSE", e);
            }
            end();
            return;
        }

        logger.info("New worker successfully started: " + this);

        boolean running = true;

        while (running && socket.isOpen() && socket.isConnected()) {
            start = System.currentTimeMillis();

            try {
                // logger.debug("wating for request...");
                request.readFrom(socket);

                // logger.debug("performing request " + request);

                if (request.getFunctionID() == AmuseMessage.FUNCTION_ID_STOP) {
                    // this will be the last call we perform
                    running = false;
                }

                if (job.hasFinished()) {
                    throw new IOException("Remote Code Proxy no longer running");
                }

                WriteMessage writeMessage = sendPort.newMessage();
                request.writeTo(writeMessage);
                writeMessage.finish();

                logger.debug("waiting for result");

                ReadMessage readMessage = null;

                while (readMessage == null) {
                    try {
                        readMessage = receivePort.receive(1000);
                    } catch (ReceiveTimedOutException timeout) {
                        // IGNORE
                    }

                    if (receivePort.connectedTo().length == 0 || job.hasFinished()) {
                        throw new IOException("receiveport no longer connected to remote proxy, or proxy no longer running");
                    }
                }
                result.readFrom(readMessage);
                readMessage.finish();

                if (result.isErrorState()) {
                    logger.warn("Error while doing call at worker", result.getError());
                }

                logger.debug("request " + request.getCallID() + " handled, result: " + result);

                // forward result to the channel
                result.writeTo(socket);

                finish = System.currentTimeMillis();

                if (logger.isDebugEnabled()) {
                    logger.debug("call took " + (finish - start) + " ms");
                }
            } catch (ConnectionClosedException e) {
                logger.info("channel closed on receiving request");
                running = false;
            } catch (IOException e) {
                running = false;
                if (socket.isOpen() && socket.isConnected()) {
                    logger.error("Error on handling call", e);

                    // report error to amuse
                    AmuseMessage errormessage =
                            new AmuseMessage(request.getCallID(), request.getFunctionID(), request.getCallCount(),
                                    "Ibis/Amuse error: " + e.getMessage());
                    try {
                        errormessage.writeTo(socket);
                    } catch (IOException e1) {
                        logger.error("Error while returning error message to amuse", e1);
                    }
                } else {
                    logger.error("Error on handling call, lost connection to AMUSE", e);
                }
            }
        }
        logger.info(this + " ending");
        end();
        logger.info(this + " done!");
    }

    private void end() {
        try {
            sendPort.close();
        } catch (IOException e) {
            logger.error("Error closing sendport", e);
        }

        try {
            job.stop();
        } catch (Exception e) {
            logger.error("Error waiting on job to finish", e);
        }

        try {
            receivePort.close(1000);
        } catch (IOException e) {
            logger.error("Error closing receiveport", e);
        }

    }

    public String toString() {
        return "Worker \"" + id;
    }
}
