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

import ibis.ipl.server.Server;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.Octopus;
import nl.esciencecenter.octopus.OctopusFactory;
import nl.esciencecenter.octopus.exceptions.OctopusException;

/**
 * Reservations of actual resources to run jobs. May still be in a queue, or already running.
 * 
 * @author Niels Drost
 * 
 */
public class ReservationManager {

    private static final Logger logger = LoggerFactory.getLogger(ReservationManager.class);

    private final ArrayList<Reservation> reservations;

    //Hubs running on frontends of resources, started if required.
    //key = resource name
    private final Map<Resource, Hub> hubs;

    private final Server iplServer;
    
    private final Octopus octopus;

    public ReservationManager(Octopus octopus) throws DistributedAmuseException {
        this.octopus = octopus;
        reservations = new ArrayList<Reservation>();
        hubs = new HashMap<Resource, Hub>();
        
        try {
            Properties properties = new Properties();
            iplServer = new Server(properties);
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create IPL server", e);
        }
        
       
    }

    public String getServerAddress() {
        return iplServer.getAddress();
    }

    public synchronized Reservation newReservation(Resource resource, String queueName, int nodeCount, int timeMinutes,
            String nodeLabel) throws DistributedAmuseException {
        logger.debug("reserving new nodes: resource name = " + resource.getName() + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        Hub hub = hubs.get(resource);

        AmuseConfiguration amuseConfiguration = resource.getConfiguration();

        if (hub == null && resource.mustStartHub()) {
            hub = new Hub(resource, amuseConfiguration, iplServer.getAddress());
            iplServer.addHubs(hub.getAddress());
            hubs.put(resource, hub);
        }

        Reservation result =
                new Reservation(resource, queueName, nodeCount, timeMinutes, nodeLabel, iplServer.getAddress(), hub,
                        amuseConfiguration);

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
