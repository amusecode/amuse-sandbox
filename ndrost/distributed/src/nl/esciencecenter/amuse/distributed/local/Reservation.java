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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.engine.util.StreamForwarder;

/**
 * @author Niels Drost
 * 
 */
public class Reservation {

    private static final Logger logger = LoggerFactory.getLogger(Reservation.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private static String buildCommand(Resource resource, AmuseConfiguration config, String serverAddress, Hub hub)
            throws DistributedAmuseException {
        StringBuilder command = new StringBuilder();

        //remote host, use ssh
        if (resource.getHostname() != null && !resource.getHostname().equals("localhost")) {
            command.append("ssh ");
            command.append(resource.getUsername());
            command.append("@");
            command.append(resource.getHostname());
            command.append(" -p ");
            command.append(resource.getPort());
            command.append(" ");
        }

        command.append(config.getJava());
        command.append(" -cp ");
        command.append(config.getAmuseHome().getPath());
        command.append("/sandbox/ndrost/distributed/distributed-server.jar");
        command.append(" nl.esciencecenter.amuse.distributed.remote.Pilot ");
        command.append(serverAddress);
        command.append(" ");
        command.append(resource.getName());

        if (hub != null) {
            command.append(" ");
            command.append(hub.getAddress());
        }

        return command.toString();
    }

    private final Process process;

    /**
     * @param resource
     * @param queueName
     * @param nodeCount
     * @param timeMinutes
     * @param nodeLabel
     * @param iplServerAddress
     * @param hub
     * @param amuseConfiguration
     */
    public Reservation(Resource resource, String queueName, int nodeCount, int timeMinutes, String nodeLabel,
            String iplServerAddress, Hub hub, AmuseConfiguration amuseConfiguration) throws DistributedAmuseException {
        try {
            String command = buildCommand(resource, amuseConfiguration, iplServerAddress, hub);

            logger.debug("running command " + command);

            process = Runtime.getRuntime().exec(command);

            new StreamForwarder(process.getErrorStream(), System.err);
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start hub on " + resource.getName(), e);
        }
    }

    public int getID() {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * 
     */
    public void cancel() {
        // TODO Auto-generated method stub

    }

    /**
     * 
     */
    public void waitUntilStarted() {
        // TODO Auto-generated method stub

    }

}
