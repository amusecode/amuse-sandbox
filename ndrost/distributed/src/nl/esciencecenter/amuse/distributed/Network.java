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

import java.util.Properties;

import ibis.ipl.Ibis;
import ibis.ipl.IbisCapabilities;
import ibis.ipl.IbisCreationFailedException;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisIdentifier;
import ibis.ipl.PortType;
import ibis.ipl.RegistryEventHandler;

/**
 * @author Niels Drost
 * 
 */
public class Network implements RegistryEventHandler {

    public static final String PORT_NAME = "amuse";

    public static PortType IPL_PORT_TYPE = new PortType(PortType.COMMUNICATION_RELIABLE, PortType.SERIALIZATION_OBJECT,
            PortType.RECEIVE_EXPLICIT, PortType.RECEIVE_TIMEOUT, PortType.CONNECTION_MANY_TO_ONE);

    public static IbisCapabilities IPL_CAPABILITIES = new IbisCapabilities(IbisCapabilities.ELECTIONS_STRICT,
            IbisCapabilities.MEMBERSHIP_TOTALLY_ORDERED, IbisCapabilities.TERMINATION, IbisCapabilities.SIGNALS);

    private final Ibis ibis;

    Network() throws DistributedAmuseException {
        try {

            Properties properties = new Properties();
            properties.put("ibis.server.address", "something only we know");
            properties.put("ibis.pool.name", "amuse");
            properties.put("ibis.location", "daemon@local");
            properties.put("ibis.managementclient", "true");
            properties.put("ibis.bytescount", "true");

            ibis = IbisFactory.createIbis(IPL_CAPABILITIES, properties, true, this, IPL_PORT_TYPE);
        } catch (IbisCreationFailedException e) {
            throw new DistributedAmuseException("failed to create ibis", e);
        }
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#died(ibis.ipl.IbisIdentifier)
     */
    @Override
    public void died(IbisIdentifier arg0) {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#electionResult(java.lang.String, ibis.ipl.IbisIdentifier)
     */
    @Override
    public void electionResult(String arg0, IbisIdentifier arg1) {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#gotSignal(java.lang.String, ibis.ipl.IbisIdentifier)
     */
    @Override
    public void gotSignal(String arg0, IbisIdentifier arg1) {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#joined(ibis.ipl.IbisIdentifier)
     */
    @Override
    public void joined(IbisIdentifier arg0) {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#left(ibis.ipl.IbisIdentifier)
     */
    @Override
    public void left(IbisIdentifier arg0) {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#poolClosed()
     */
    @Override
    public void poolClosed() {
        // TODO Auto-generated method stub
        
    }

    /* (non-Javadoc)
     * @see ibis.ipl.RegistryEventHandler#poolTerminated(ibis.ipl.IbisIdentifier)
     */
    @Override
    public void poolTerminated(IbisIdentifier arg0) {
        // TODO Auto-generated method stub
        
    }

    /**
     * @return
     */
    public Ibis getIbis() {
        // TODO Auto-generated method stub
        return null;
    }

}
