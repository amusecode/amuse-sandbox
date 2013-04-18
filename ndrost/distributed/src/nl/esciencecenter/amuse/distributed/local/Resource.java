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

/**
 * @author Niels Drost
 * 
 */
public class Resource {

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

    public Resource(String name, String hostname, String amuseDir, int port, String username, String schedulerType) {
        this.id = getNextID();
        this.name = name;
        this.hostname = hostname;
        this.amuseDir = amuseDir;
        this.port = port;
        this.username = username;
        this.schedulerType = schedulerType;
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
    
}
