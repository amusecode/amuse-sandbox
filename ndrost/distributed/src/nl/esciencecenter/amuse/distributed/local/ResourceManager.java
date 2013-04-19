package nl.esciencecenter.amuse.distributed.local;

import java.util.ArrayList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;

/**
 * Manages resources potentially available for starting reservations on.
 * 
 * @author Niels Drost
 * 
 */
public class ResourceManager {

    private static final Logger logger = LoggerFactory.getLogger(ResourceManager.class);

    private ArrayList<Resource> resources;

    ResourceManager() {
        resources = new ArrayList<Resource>();
    }

    public synchronized Resource newResource(String name, String hostname, String amuseDir, int port, String username,
            String schedulerType, Boolean startHub) throws DistributedAmuseException {
        logger.debug("creating new resource: name = " + name + " hostname = " + hostname + " port = " + port + " user name = "
                + username + " scheduler type = " + schedulerType + " amuse dir = " + amuseDir + " start hub = " + startHub);

        for (Resource resource : resources) {
            if (resource.getName().equals(name)) {
                throw new DistributedAmuseException("Resource " + name + " already exists");
            }
        }

        Resource result = new Resource(name, hostname, amuseDir, port, username, schedulerType, startHub);

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
}
