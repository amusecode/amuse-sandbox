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


import java.util.ArrayList;
import java.util.Iterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.IbisIdentifier;
import ibis.ipl.RegistryEventHandler;

/**
 * Set of nodes available to run jobs on
 * 
 * @author Niels Drost
 *
 */
public class Nodes implements RegistryEventHandler {
    
    private static final Logger logger = LoggerFactory.getLogger(Nodes.class);
    
    private final ArrayList<Node> nodes;
    
    /**
     * @param jobManager
     */
    public Nodes() {
        nodes = new ArrayList<Node>();
    }
    
    /**
     * @return
     */
    public synchronized boolean isEmpty() {
        return nodes.isEmpty();
    }

    /**
     * @param job
     * @return
     */
    public synchronized Node[] getSuitableNodes(Job job) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public synchronized void died(IbisIdentifier ibis) {
        //handle like it left
        left(ibis);
    }

    @Override
    public void joined(IbisIdentifier ibis) {
        logger.debug("new Ibis joined: " + ibis);
        
        nodes.add(new Node(ibis));
    }

    @Override
    public void left(IbisIdentifier ibis) {
        logger.debug("Ibis left: " + ibis);
        
        Iterator<Node> iterator = nodes.iterator();
        
        while(iterator.hasNext()) {
            if (iterator.next().getIbisIdentifier().equals(ibis)) {
                iterator.remove();
                //TODO: do something with the jobs still running on this node, if any...
            }
        }
    }
    
    @Override
    public void electionResult(String name, IbisIdentifier winner) {
        //IGNORED
    }

    @Override
    public void gotSignal(String signal, IbisIdentifier origin) {
        //IGNORED
    }

    @Override
    public void poolClosed() {
        //IGNORED
    }

    @Override
    public void poolTerminated(IbisIdentifier arg0) {
        //IGNORED
    }

  

}
