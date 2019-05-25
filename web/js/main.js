var cytoscape = require('cytoscape');
var dagre = require('cytoscape-dagre');
var data = require('./networks');
var path = require('path');
var express = require('express');
cytoscape.use(dagre);

const app = express();
const port = 3000

var cy = cytoscape({
    elements: data,
    // interaction options:
    minZoom: 1e-50,
    maxZoom: 1e50,
    zoomingEnabled: true,
    userZoomingEnabled: true,
    panningEnabled: true,
    userPanningEnabled: true,
    boxSelectionEnabled: false,
    selectionType: 'single',
    touchTapThreshold: 8,
    desktopTapThreshold: 4,
    autolock: false,
    autoungrabify: false,
    autounselectify: false,

    headless: true,
});


cy.layout({
    name: 'dagre',
    defaults : {
        // dagre algo options, uses default value on undefined
        nodeSep: undefined, // the separation between adjacent nodes in the same rank
        edgeSep: undefined, // the separation between adjacent edges in the same rank
        rankSep: undefined, // the separation between adjacent nodes in the same rank
        rankDir: undefined, // 'TB' for top to bottom flow, 'LR' for left to right,
        ranker: undefined, // Type of algorithm to assign a rank to each node in the input graph. Possible values: 'network-simplex', 'tight-tree' or 'longest-path'
        minLen: function( edge ){ return 1; }, // number of ranks to keep between the source and target of the edge
        edgeWeight: function( edge ){ return 1; }, // higher weight edges are generally made shorter and straighter than lower weight edges
      
        // general layout options
        fit: true, // whether to fit to viewport
        padding: 30, // fit padding
        spacingFactor: undefined, // Applies a multiplicative factor (>0) to expand or compress the overall area that the nodes take up
        nodeDimensionsIncludeLabels: false, // whether labels should be included in determining the space used by a node
        animate: false, // whether to transition the node positions
        animateFilter: function( node, i ){ return true; }, // whether to animate specific nodes when animation is on; non-animated nodes immediately go to their final positions
        animationDuration: 500, // duration of animation in ms if enabled
        animationEasing: undefined, // easing of animation if enabled
        boundingBox: undefined, // constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
        transform: function( node, pos ){ return pos; }, // a function that applies a transform to the final node position
        ready: function(){}, // on layoutready
        stop: function(){} // on layoutstop
      },
}).run();


app.get('/', (req, res) => res.sendFile(path.resolve("../index.html")));
//app.get('/', (req, res) => res.send("hello"));

app.listen(port, () => console.log(`app listening on port ${port}!`));