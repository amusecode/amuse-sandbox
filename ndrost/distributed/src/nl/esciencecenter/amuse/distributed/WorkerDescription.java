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
 * Description of a worker
 * 
 * @author Niels Drost
 */
public class WorkerDescription {

    private final String id;
    private final String codeName;
    private final String codeDir;
    private final String stdoutFile;
    private final String stderrFile;
    private final String nodeLabel;

    private final int nrOfWorkers;
    private final int nrOfNodes;
    private final int nrOfThreads;

    private final boolean copyWorkerCode;

    public WorkerDescription(String id, String codeName, String codeDir, String stdoutFile, String stderrFile, String nodeLabel,
            int nrOfWorkers, int nrOfNodes, int nrOfThreads, boolean copyWorkerCode) {
        this.id = id;
        this.codeName = codeName;
        this.codeDir = codeDir;
        this.stdoutFile = stdoutFile;
        this.stderrFile = stderrFile;
        this.nodeLabel = nodeLabel;
        this.nrOfWorkers = nrOfWorkers;
        this.nrOfNodes = nrOfNodes;
        this.nrOfThreads = nrOfThreads;
        this.copyWorkerCode = copyWorkerCode;
    }

    public String getID() {
        return id;
    }

    public String getCodeName() {
        return codeName;
    }

    public String getCodeDir() {
        return codeDir;
    }

    public String getStdoutFile() {
        return stdoutFile;
    }

    public String getStderrFile() {
        return stderrFile;
    }

    public String getNodeLabel() {
        return nodeLabel;
    }

    public int getNrOfWorkers() {
        return nrOfWorkers;
    }

    public int getNrOfNodes() {
        return nrOfNodes;
    }

    public int getNrOfThreads() {
        return nrOfThreads;
    }

    public boolean copyWorkerCode() {
        return copyWorkerCode;
    }
}
