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
package nl.esciencecenter.amuse.distributed;

/**
 * @author Niels Drost
 *
 */
public class Reservation {
    
    private static int nextID = 0;

    /**
     * @param resource
     * @param queueName
     * @param nodeCount
     * @param timeMinutes
     * @param nodeLabel
     */
    public Reservation(Resource resource, String queueName, int nodeCount, int timeMinutes, String nodeLabel) {
        // TODO Auto-generated constructor stub
    }

    private static int getNextID() {
        return nextID++;
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
