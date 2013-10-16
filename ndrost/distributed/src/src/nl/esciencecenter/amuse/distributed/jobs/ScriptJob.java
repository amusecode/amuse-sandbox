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

import ibis.ipl.Ibis;
import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.IOException;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

/**
 * @author Niels Drost
 *
 */
public class ScriptJob extends Job {

    /**
     * @param nodeLabel
     * @param numberOfNodes
     * @param ibis
     * @throws DistributedAmuseException 
     * @throws IOException 
     */
    public ScriptJob(String nodeLabel, int numberOfNodes, Ibis ibis) throws DistributedAmuseException {
        super(nodeLabel, numberOfNodes, ibis);
        // TODO Auto-generated constructor stub
    }

    /**
     * @param writeMessage
     * @throws IOException
     */
    @Override
    void writeJobDetails(WriteMessage writeMessage) throws IOException {
        //FIXME: transfer files etc
        // TODO Auto-generated method stub
    }

    /**
     * @param readMessage
     * @throws ClassNotFoundException
     * @throws IOException
     */
    @Override
    void readJobStatus(ReadMessage readMessage) throws ClassNotFoundException, IOException {
        // TODO Auto-generated method stub
        
    }

    /**
     * @param readMessage
     * @throws ClassNotFoundException
     * @throws IOException
     */
    @Override
    void readJobResult(ReadMessage readMessage) throws ClassNotFoundException, IOException {
        // TODO Auto-generated method stub
        
    }

    /**
     * @return
     */
    public boolean useCodeCache() {
        // TODO Auto-generated method stub
        return false;
    }

    /**
     * @return
     */
    public String getScriptName() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * @return
     */
    public String getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * @return
     */
    public String getScriptDir() {
        // TODO Auto-generated method stub
        return null;
    }

}
