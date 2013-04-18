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

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.Network;

import java.util.ArrayList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DistributedAmuse {

    private static final Logger logger = LoggerFactory.getLogger(DistributedAmuse.class);

    //resources potentially available for starting reservations on
    private final ArrayList<Resource> resources;

    //reservations of actual resources to run jobs. May still be in a queue, or already running
    private final ArrayList<Reservation> reservations;

    //takes care of job queue
    private final AmuseJobScheduler scheduler;

    private final WorkerConnectionHandler workerConnectionHandler;

    private final Network network;

    public DistributedAmuse() throws DistributedAmuseException {
        resources = new ArrayList<Resource>();
        reservations = new ArrayList<Reservation>();

        network = new Network();

        scheduler = new AmuseJobScheduler(this);

        workerConnectionHandler = new WorkerConnectionHandler(this);
    }

    public Network getNetwork() {
        return network;
    }
    
    public AmuseJobScheduler getScheduler() {
        return scheduler;
    }

    /**
     * Port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     * 
     * @return the port used by the IbisChannel to connect to when creating workers and stderr/stdout streams
     */
    public synchronized int getWorkerPort() {
        logger.debug("returning worker port");
        return workerConnectionHandler.getPort();
    }

    public synchronized Resource newResource(String name, String hostname, String amuseDir, int port, String username,
            String schedulerType) throws DistributedAmuseException {
        logger.debug("creating new resource: name = " + name + " hostname = " + hostname + " port = " + port + " user name = "
                + username + " scheduler type = " + schedulerType + " amuse dir = " + amuseDir);

        for (Resource resource : resources) {
            if (resource.getName().equals(name)) {
                throw new DistributedAmuseException("Resource " + name + " already exists");
            }
        }

        Resource result = new Resource(name, hostname, amuseDir, port, username, schedulerType);

        resources.add(result);

        return result;
    }

    public synchronized Resource getResource(int resourceID) throws DistributedAmuseException {
        for (Resource resource : resources) {
            if (resource.getId() == resourceID) {
                return resource;
            }
        }
        throw new DistributedAmuseException("Resource with ID " + resourceID + " not found");
    }

    public synchronized Resource getResource(String name) throws DistributedAmuseException {
        for (Resource resource : resources) {
            if (resource.getName().equals(name)) {
                return resource;
            }
        }
        throw new DistributedAmuseException("Resource with name " + name + " not found");
    }

    public synchronized void deleteResource(Resource resource) throws DistributedAmuseException {
        for (int i = 0; i < resources.size(); i++) {
            if (resource.getId() == resources.get(i).getId()) {
                resources.remove(i);
                return;
            }
        }
        throw new DistributedAmuseException("Resource " + resource.getId() + " not found");
    }

    public synchronized Reservation newReservation(String resourceName, String queueName, int nodeCount, int timeMinutes,
            String nodeLabel) throws DistributedAmuseException {
        logger.debug("reserving new nodes: resource name = " + resourceName + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        Resource resource = getResource(resourceName);

        Reservation result = new Reservation(resource, queueName, nodeCount, timeMinutes, nodeLabel);

        reservations.add(result);

        return result;
    }

    public synchronized Reservation getReservation(int reservationID) throws DistributedAmuseException {
        for (Reservation reservation : reservations) {
            if (reservation.getID() == reservationID) {
                return reservation;
            }
        }
        throw new DistributedAmuseException("Reservation with ID " + reservationID + " not found");
    }

    public synchronized void deleteReservation(int reservationID) throws DistributedAmuseException {
        logger.debug("deleting reservation " + reservationID);

        for (int i = 0; i < reservations.size(); i++) {
            Reservation reservation = reservations.get(i);
            if (reservationID == reservation.getID()) {
                reservations.remove(i);
                reservation.cancel();
                return;
            }
        }
        throw new DistributedAmuseException("Reservation " + reservationID + " not found");
    }

    public synchronized void waitForAllReservations() {
        logger.debug("waiting for all reservations to start");

        for (Reservation reservation : reservations) {
            reservation.waitUntilStarted();
        }
    }

}
