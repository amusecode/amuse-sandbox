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
package nl.esciencecenter.amuse.distributed.remote;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;

import java.util.Properties;

import nl.esciencecenter.amuse.distributed.Network;

/**
 * Pilot job. Started when a reservations is made.
 * 
 * @author Niels Drost
 * 
 */
public class Pilot {

    public static void main(String[] arguments) throws Throwable {
        Properties properties = new Properties();
        properties.put("ibis.location", arguments[0]);
        properties.put("ibis.server.address", arguments[1]);
        properties.put("ibis.pool.name", "amuse");

        if (arguments.length == 3) {
            properties.put("ibis.hub.addresses", arguments[2]);
        }
        //properties.put("ibis.managementclient", "true");
        //properties.put("ibis.bytescount", "true");

        Ibis ibis = IbisFactory.createIbis(Network.IPL_CAPABILITIES, properties, true, null, Network.IPL_PORT_TYPE);

        System.err.println("running Pilot at location \"" + arguments[0] + "\" for 60 seconds");
        
        Thread.sleep(60000);
        
        ibis.end();
    }
}
