<script src="//d3js.org/d3.v3.min.js" charset="utf-8"></script>

var force = d3.layout.force()
  .charge(-3000)
  .linkDistance(50)
  .theta(0.2)
  .gravity(1)
  .size([width, height]);


function GO(goid) {
  return goNodes[goid];
}


 //drawing the graph
function drawGraph(struct) {
  force
    .nodes(struct.nodes)
    .links(struct.links);

  var strokeColor = d3.scale.category10();
  var radiusScale = d3.scale.linear().domain([0,d3.max(struct.nodes,function(d) {return d.r;})]).range([5,10]);

  var link = svg.selectAll("line.link")
    .data(struct.links)
    .enter().append("line")
    .attr("class", "link")
    .style('stroke', function(d, i ) { return strokeColor(d.type);})
    .style("stroke-width", 2);
 
  var node = svg.selectAll("circle.node")
    .data(struct.nodes)
    .enter().append("circle")
    .attr("class", "node")
    .attr("r", function(d) { return radiusScale(d.r);})
    .style("fill", function(d,i) {
      if (d.is == "root") {
        return 'black';
      } else if (d.is == "target" || d.is == 'child') {
        return 'orange';
      } else {
        return 'lightblue';
      }
    })
    .call(force.drag)
    .on('click',function(d,i) {
      //update(d['id']);
    });
 
  var nodeText = svg.selectAll("text.node")
    .data(struct.nodes)
    .enter().append("text")
    .attr("class", "node")
    .attr("displayState",'y')
    //.style('display','none')
    .text(function(d) { return d.name;});
 
  force.on("tick", function(e) {
    link.attr("x1", function(d) { return d.source.x; })
      .attr("y1", function(d) { return d.source.y; })
      .attr("x2", function(d) { return d.target.x; })
      .attr("y2", function(d) { return d.target.y; });
 
    node.attr("cx", function(d) { return d.x; })
      .attr("cy", function(d) { return d.y; })

    nodeText.attr("x", function(d) { return d.x; })
      .attr("y", function(d) { return d.y + radiusScale(d.r) + 10; })
  });

  force.start();
}


function getUniqueGoid(goids){
  var dic={},r=[];
  for(i=0;i<goids.length;i++){
    console.log(goids[i][0]);
    dic[goids[i][0]] = goids[i][1];
  }
  for(i in dic){
    arr = [i,dic[i]];
    r.push(arr);
  }
  return r;
}



//creates a D3 structure that can be inputted into the force-directed graph layout
function createD3Structure(goids, target, width , height) {

  goids = getUniqueGoid(goids);

  nodes = [];
  links = [];
  
  nodeIndex = {};
  levels = goids[goids.length-1][1]; //total level
  level_width ={};

  for (var i = 0, count = goids.length ; i < count ; i ++ ) {
    go_level = goids[i];
    goid = go_level[0];
    level = go_level[1];
    node = {};
    node['name'] = GO(goid)['n'];
    node['id'] = goid;
    node['r'] = 10;
    node['is'] = 'node';
    nodeIndex[goid] = i;
    

    //get the number of GO on each level
    if(level_width[level]==undefined){
      level_width[level]=1;
      node['position'] = 1;
    }
    else{
      level_width[level]++;
      node['position'] = level_width[level];
    }

    nodes.push(node);
  }

  for (var i = 0, count = goids.length ; i < count ; i ++ ) {
    go_level = goids[i];
    goid = go_level[0];
    level = go_level[1];

    parents = GO(goid)['p'];

    if (parents.length > 0 ) {
      for (var a = 0, countA = parents.length ; a < countA ; a ++ ) {
        parent = parents[a];
        pgoid = parent;
        pIndex = nodeIndex[pgoid];

     /*   position = nodes[pIndex]['position']; 
        nodes[pIndex].fixed = true;
        nodes[pIndex].x = width*(position/level_width[level]*0.5);
        nodes[pIndex].y = height*(1-(level)/levels);*/

        if (pgoid == target) {
          nodes[i]['is'] = 'child';
          nodes[i]['r'] = 30;
        }

        if (pIndex != undefined) {
          link = {};
          link['source'] = pIndex;
          link['target'] = i;
          link['value'] = 5;
          //link['type'] = prel;
          links.push(link)
        }
      }
    } else {        
      nodes[i].fixed = true; //the top node
      nodes[i].x = width / 2;
      nodes[i].y = 50;
      nodes[i]['is'] = 'root';
    }
  }


  nodes[nodeIndex[target]].fixed = true;
  nodes[nodeIndex[target]].x = width / 2;
  nodes[nodeIndex[target]].y = height - 20;
  nodes[nodeIndex[target]]['is'] = 'target';

  return {'nodes' : nodes, 'links' : links};
}      