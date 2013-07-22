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
package nl.esciencecenter.amuse.distributed.jobs;

import ibis.ipl.IbisIdentifier;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Node running one or more jobs. Note: NOT Thread-safe.
 * 
 * @author Niels Drost
 * 
 */
public class Node {

    private static final Logger logger = LoggerFactory.getLogger(Node.class);

    private int slots;

    private final String label;

    //address of this node
    private final IbisIdentifier ibisIdentifier;

    //list of all jobs running on this node
    private List<Job> jobs;

    public Node(IbisIdentifier ibisIdentifier) {
        this.ibisIdentifier = ibisIdentifier;

        jobs = new LinkedList<Job>();

        String[] tags = ibisIdentifier.tagAsString().split(",");

        if (tags.length != 2) {
            logger.error("Cannot parse ibis tag: " + ibisIdentifier.tagAsString());
            label = "unknown";
            slots = 1;
        } else {

            label = tags[0];

            try {
                slots = Integer.parseInt(tags[1]);
            } catch (NumberFormatException e) {
                logger.error("Cannot parse ibis tag: " + ibisIdentifier.tagAsString(), e);
                slots = 1;
            }
        }

    }

    void addJob(Job job) {
        jobs.add(job);
    }

    boolean isAvailableForBatchJobs() {
        int batchJobCount = 0;
        int workerJobCount = 0;

        Iterator<Job> iterator = jobs.iterator();

        while (iterator.hasNext()) {
            Job job = iterator.next();

            if (job.isDone()) {
                iterator.remove();
            } else if (job.isBatchJob()) {
                batchJobCount++;
            } else if (job.isWorkerJob()) {
                workerJobCount++;
            }
        }

        return (batchJobCount < slots) && (workerJobCount == 0);
    }

    public int getSlots() {
        return slots;
    }

    public String getLabel() {
        return label;
    }

    public IbisIdentifier getIbisIdentifier() {
        return ibisIdentifier;
    }
}
