package ibis.amuse;

import java.io.Serializable;

/**
 * Settings for a worker.
 * 
 * @author Niels Drost
 * 
 */
public class WorkerInfo implements Serializable {

    private static final long serialVersionUID = 1L;

    private String id;

    private int nrOfNodes;

    private int nrOfProcesses;

    private String codeName;

    private String codeDir;

    private String amuseHome;

    private String mpiexec;

    private String mpdboot;

    private String interfaceType;

    String getID() {
        return id;
    }

    void setWorkerID(String workerID) {
        this.id = workerID;
    }

    int getNrOfNodes() {
        return nrOfNodes;
    }

    void setNrOfNodes(int nrOfNodes) {
        this.nrOfNodes = nrOfNodes;
    }

    int getNrOfProcesses() {
        return nrOfProcesses;
    }

    void setNrOfProcesses(int nrOfProcesses) {
        this.nrOfProcesses = nrOfProcesses;
    }

    String getCodeName() {
        return codeName;
    }

    void setCodeName(String codeName) {
        this.codeName = codeName;
    }

    String getCodeDir() {
        return codeDir;
    }

    void setCodeDir(String codeDir) {
        this.codeDir = codeDir;
    }

    String getAmuseHome() {
        return amuseHome;
    }

    void setAmuseHome(String amuseHome) {
        this.amuseHome = amuseHome;
    }

    String getMpiexec() {
        return mpiexec;
    }

    void setMpiexec(String mpiexec) {
        this.mpiexec = mpiexec;
    }

    String getMpdboot() {
        return mpdboot;
    }

    void setMpdboot(String mpdboot) {
        this.mpdboot = mpdboot;
    }

    String getInterfaceType() {
        return interfaceType;
    }

    void setInterfaceType(String interfaceType) {
        this.interfaceType = interfaceType;
    }

}
