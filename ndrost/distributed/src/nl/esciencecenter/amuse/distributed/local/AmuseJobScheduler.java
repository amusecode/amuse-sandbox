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
import nl.esciencecenter.amuse.distributed.WorkerDescription;

/**
 * @author Niels Drost
 * 
 */
public class AmuseJobScheduler {

    /**
     * @param distributedAmuse
     */
    public AmuseJobScheduler(DistributedAmuse distributedAmuse) {
        // TODO Auto-generated constructor stub
    }

    /**
     * @param job
     */
    public void waitForJob(int jobID) throws DistributedAmuseException {
        // TODO Auto-generated method stub

    }

    /**
     * @param jobID
     * @return
     */
    public String getJobResult(int jobID) throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return "result!";
    }

    /**
     * 
     */
    public void waitForAllJobs() throws DistributedAmuseException {
        // TODO Auto-generated method stub

    }

    /**
     * @param script
     * @param arguments
     * @param codeDir
     * @param nodeLabel
     * @param useCodeCache
     * @return
     */
    public int submitScriptJob(String script, String arguments, String codeDir, String nodeLabel, boolean useCodeCache)
            throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * @param function
     * @param arguments
     * @param nodeLabel
     * @return
     */
    public int submitPickledJob(String function, String arguments, String nodeLabel) throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return 0;
    }

    /**
     * @param workerDescription
     * @return
     */
    public AmuseWorkerJob submitWorkerJob(WorkerDescription workerDescription) throws DistributedAmuseException {
        // TODO Auto-generated method stub
        return null;
    }

}
