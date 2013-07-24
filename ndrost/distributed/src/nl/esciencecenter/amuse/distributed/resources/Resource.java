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
package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.Server;

import java.io.InputStream;
import java.net.URI;
import java.util.Arrays;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.Octopus;
import nl.esciencecenter.octopus.credentials.Credential;
import nl.esciencecenter.octopus.files.AbsolutePath;
import nl.esciencecenter.octopus.files.FileSystem;
import nl.esciencecenter.octopus.files.RelativePath;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 */
public class Resource {

    private static final Logger logger = LoggerFactory.getLogger(Resource.class);

    private static int nextID = 0;

    private static int getNextID() {
        return nextID++;
    }

    private final int id;
    private final String name;
    private final String hostname;
    private final String amuseDir;
    private final int port;
    private final String username;
    private final String schedulerType;

    private final AmuseConfiguration configuration;

    //null == auto
    private final Boolean startHub;

    private final Hub hub;
    
    private static void waitUntilHubStarted(Server iplServer, String hubAddress, String name) throws DistributedAmuseException {
        for (int i = 0; i < 40; i++) {
            String[] knownHubAddresses = iplServer.getHubs();
            logger.debug("ipl hub addresses now " + Arrays.toString(iplServer.getHubs()));
            for (String knownHub : knownHubAddresses) {
                if (knownHub.equals(hubAddress)) {
                    logger.debug("new hub now connected to local hub");
                    return;
                }
            }
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                //IGNORE
            }
        }
        throw new DistributedAmuseException("Local and new remote Hub at " + name + " not able to communicate");
    }

    public Resource(String name, String hostname, String amuseDir, int port, String username, String schedulerType,
            Boolean startHub, Octopus octopus, Server iplServer) throws DistributedAmuseException {
        this.id = getNextID();
        this.name = name;
        this.hostname = hostname;
        this.amuseDir = amuseDir;
        this.port = port;
        this.username = username;
        this.schedulerType = schedulerType;
        this.startHub = startHub;

        this.configuration = downloadConfiguration(octopus);

        if (mustStartHub()) {
            this.hub = new Hub(this, this.configuration, iplServer.getHubs(), octopus);
            iplServer.addHubs(this.hub.getAddress());

            String hubAddress = this.hub.getAddress();
            logger.debug("just added new hub " + hubAddress);

            waitUntilHubStarted(iplServer, hubAddress, name);
        } else {
            this.hub = null;
        }
        logger.debug("Created new resource {}", this);
    }

    private AmuseConfiguration downloadConfiguration(Octopus octopus) throws DistributedAmuseException {
        try {
            FileSystem filesystem;

            if (this.name.equals("local")) {
                filesystem = octopus.files().getLocalHomeFileSystem();
            } else {
                URI uri = new URI("ssh", username, hostname, port, null, null, null);
                Credential credential = octopus.credentials().getDefaultCredential("ssh");

                filesystem = octopus.files().newFileSystem(uri, credential, null);
            }

            RelativePath amuseConfig = new RelativePath(this.amuseDir + "/config.mk");

            AbsolutePath path = octopus.files().newPath(filesystem, amuseConfig);

            try (InputStream in = octopus.files().newInputStream(path)) {
                return new AmuseConfiguration(this.amuseDir, in);
            }
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot download configuration file for resource " + this.name, e);
        }
    }

    public int getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public String getHostname() {
        return hostname;
    }

    public String getAmuseDir() {
        return amuseDir;
    }

    public int getPort() {
        return port;
    }

    public String getUsername() {
        return username;
    }

    public String getSchedulerType() {
        return schedulerType;
    }

    public AmuseConfiguration getConfiguration() {
        return configuration;
    }

    @Override
    public int hashCode() {
        return new Integer(id).hashCode();
    }

    @Override
    public boolean equals(Object other) {
        if (other == null) {
            return false;
        }

        if (!(other instanceof Resource)) {
            return false;
        }

        return id == ((Resource) other).id;
    }

    public boolean mustStartHub() {
        if (startHub == null) {
            return getSchedulerType() != null && (getSchedulerType().toLowerCase() != "local");
        }
        return startHub;
    }

    public void stop() {
        logger.debug("Stopping resource {}", this);
        if (hub != null) {
            hub.stop();
        }
    }

    public Hub getHub() {
        return hub;
    }

    public boolean isLocal() {
        return hostname == null || hostname.equals("localhost") || hostname.equals("local");
    }

    @Override
    public String toString() {
        return "Resource [id=" + id + ", name=" + name + ", hostname=" + hostname + ", amuseDir=" + amuseDir + ", port=" + port
                + ", username=" + username + ", schedulerType=" + schedulerType + ", configuration=" + configuration
                + ", startHub=" + startHub + ", hub=" + hub + "]";
    }
}
