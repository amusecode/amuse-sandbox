/**
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the 'License'); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an 'AS IS' BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
$(document).ready(function(){

  var client, destination;

    var url = 'ws://localhost:61623/stomp';
    var login = 'guest';
    var passcode = 'password';
    destination = '/queue/data';
    
    var plot = $.plot(
        "#messages", 
        [ 
            {
                data: [[0,0]],
                points: { show: true }
            } 
        ], 
        {
                series: {
                    shadowSize: 0	// Drawing is faster without shadows
                },
                yaxis: {
                    min: -2,
                    max: 2
                },
                xaxis: {
                    min: -2,
                    max: 2
                }
        }
    );

    var queue = [];
    setInterval(function() {
        if(queue.length == 0) {
            return;
        }
        var message = queue.shift();
        
        if(message.type == 'particles') {
            var points = [];
            for (var i = 0; i < message.x.length; ++i) {
                points.push([message.x[i],message.y[i]]);
            }
            plot.setData([{data:points, points:{show:true}}]);
            plot.draw()
            $('#title').text(message['time-str']);
        } else if (message.type == 'encounter'){
            var table = $('#encounters');
            var row = $('<tr></tr>');
            var before = $('<td></td>');
            before.addClass('encounter-view');
            var after = $('<td></td>');
            after.addClass('encounter-view');
            row.append(before);
            row.append(after);
            table.append(row);
            var before_data = [];
            for (var i = 0; i < message.before.length; ++i) {
                before_data.push([message.before[i].x,message.before[i].y,message.before[i].vx,message.before[i].vy]);
            }
            var p = $.plot(
                before, 
                [ 
                    {
                        data: before_data,
                        points: { show: true }
                    } 
                ], 
                {
                    series: {
                        shadowSize: 0	// Drawing is faster without shadows
                        ,
                        direction: {
                            show: true
                        }
                    }
                }   
            );
            var after_data = [];
            for (var i = 0; i < message.after.length; ++i) {
                after_data.push([message.after[i].x,message.after[i].y,message.after[i].vx,message.after[i].vy]);
            }
            var p2 = $.plot(
                after, 
                [ 
                    {
                        data: after_data,
                        points: { show: true }
                    } 
                ], 
                {
                     series: {
                        shadowSize: 0	// Drawing is faster without shadows
                        ,
                        direction: {
                            show: true
                        }
                    }
                }   
            );
            
        }
    }, 50);
    
    client = Stomp.client(url); //, ['binary']);

    // this allows to display debug logs directly on the web page
    client.debug = function(str) {
      //$("#debug").append(document.createTextNode(str + "\n"));
      //console.log(str);
    };
    // the client is notified when it is connected to the server.
    var onconnect = function(frame) {
        client.debug("connected to Stomp");
        client.subscribe(
            destination, 
            function(message) {
                var message = JSON.parse(message.body);
                queue.push(message);
            }
            ,
            {
                'browser':true,
                'browser-end':false,
                'include-seq':'seq',
                'from-seq':0
            }
        );
    };
    client.connect(login, passcode, onconnect);


});
