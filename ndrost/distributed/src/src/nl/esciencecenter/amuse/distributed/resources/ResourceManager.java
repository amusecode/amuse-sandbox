package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.Server;
import ibis.ipl.server.ServerProperties;

import java.io.File;
import java.util.ArrayList;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.Octopus;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Manages resources potentially available for starting reservations on.
 * 
 * @author Niels Drost
 * 
 */
public class ResourceManager {

    private static final Logger logger = LoggerFactory.getLogger(ResourceManager.class);

    private final Octopus octopus;
    
    private final Server iplServer;
    
    private final ArrayList<Resource> resources;
    
    public ResourceManager(Octopus octopus, File tmpDir, String amuseRootDir) throws DistributedAmuseException {
        resources = new ArrayList<Resource>();
        this.octopus = octopus;
        
        try {
            Properties properties = new Properties();
            //use a random free port.
            properties.put(ServerProperties.PORT, "0");
            this.iplServer = new Server(properties);
        } catch (Exception e) {
            throw new DistributedAmuseException("could not create IPL server", e);
        }

        //add local resource by default

        logger.debug("local amuse dir = " + amuseRootDir);
        newResource("local", null, amuseRootDir, -1, null, "local", false);
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

        Resource result = new Resource(name, hostname, amuseDir, port, username, schedulerType, startHub, octopus, iplServer);

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
    
    public synchronized int getResourceCount() {
        return resources.size();
    }
    
    public synchronized Resource[] getResources() {
        return resources.toArray(new Resource[resources.size()]);
    }

    public synchronized void deleteResource(Resource resource) throws DistributedAmuseException {
        for (int i = 0; i < resources.size(); i++) {
            if (resource.getId() == resources.get(i).getId()) {
                resource.stop();
                resources.remove(i);
                return;
            }
        }
        throw new DistributedAmuseException("Resource " + resource.getId() + " not found");
    }

    public String getIplServerAddress() {
        return iplServer.getAddress();
    }

    public String[] getHubAddresses() {
        return iplServer.getHubs();
    }
    
    public synchronized void end() {
        for(Resource resource: resources) {
            resource.stop();
        }
        iplServer.end(1000);
    }
}
