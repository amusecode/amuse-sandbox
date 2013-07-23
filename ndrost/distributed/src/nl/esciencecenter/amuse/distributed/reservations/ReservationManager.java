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
package nl.esciencecenter.amuse.distributed.reservations;

import java.util.ArrayList;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.resources.Resource;
import nl.esciencecenter.amuse.distributed.resources.ResourceManager;
import nl.esciencecenter.octopus.Octopus;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Reservations of actual resources to run jobs. May still be in a queue, or already running.
 * 
 * @author Niels Drost
 * 
 */
public class ReservationManager {

    private static final Logger logger = LoggerFactory.getLogger(ReservationManager.class);

    private final ResourceManager resourceManager;

    private final ArrayList<Reservation> reservations;

    private final Octopus octopus;

    public ReservationManager(Octopus octopus, ResourceManager resourceManager) throws DistributedAmuseException {
        this.octopus = octopus;
        this.resourceManager = resourceManager;
        reservations = new ArrayList<Reservation>();
    }

    public synchronized Reservation newReservation(String resourceName, String queueName, int nodeCount, int timeMinutes,
            String nodeLabel) throws DistributedAmuseException {
        logger.debug("reserving new nodes: resource name = " + resourceName + " queue name = " + queueName
                + " number of nodes = " + nodeCount + " time (in minutes) = " + timeMinutes + " node label = " + nodeLabel);

        Resource resource = resourceManager.getResource(resourceName);

        Reservation result =
                new Reservation(resource, queueName, nodeCount, timeMinutes, nodeLabel, resourceManager.getIplServerAddress(),
                        resourceManager.getHubAddresses(), octopus);

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

    public synchronized void waitForAllReservations() throws DistributedAmuseException {
        logger.debug("waiting for all reservations to start");

        for (Reservation reservation : reservations) {
            reservation.waitUntilStarted();
        }

        logger.debug("All reservations started");
    }

    public synchronized void end() {
        for (Reservation reservation : reservations) {
            try {
                reservation.cancel();
            } catch (DistributedAmuseException e) {
                logger.error("Failed to cancel reservation: " + reservation, e);
            }
        }
    }
}
