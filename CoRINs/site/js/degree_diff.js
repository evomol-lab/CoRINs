var pathAuxDegreeDiff = "";

function degree_diff(path , widthGraf = window.innerWidth * 0.40 ){

    pathAuxDegreeDiff = path;

    var margin = {top: 40, right: 20, bottom: 100, left: 40},
        //width = 960 - margin.left - margin.right,
        width = widthGraf - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;

    var formatPercent = d3v3.format("1");

    var x = d3v3.scale.ordinal()
        .rangeRoundBands([0, width], .1);

    var y = d3v3.scale.linear()
        .range([height, 0]);

    var xAxis = d3v3.svg.axis()
        .scale(x)
        .orient("bottom");    

    var yAxis = d3v3.svg.axis()
        .scale(y)
        .orient("left")
        .tickFormat(formatPercent);

    /*var tip = d3v3.tip()
    .attr('class', 'd3-tip')
    .attr('id', 'd3tip')
    .offset([-10, 0])
    .html(function(d) {
        return "<strong>Differences:</strong> <span style='color:red'>" + d.differences + "</span>";
    })*/

  var Tooltip = d3v3.select("#degree_diff")
  .append("div")
  .style("opacity", 0)
  .attr("class", "tooltip")
  .style("background-color", "white")
  .style("border", "solid")
  .style("border-width", "2px")
  .style("border-radius", "5px")
  .style("padding", "5px")

// Three function that change the tooltip when user hover / move / leave a cell
var mouseover = function(d) {
  Tooltip
    .style("opacity", 1)
  d3v3.select(this)
    .style("stroke", "black")
    .style("opacity", 1)
}
var mousemove = function(d) {
  Tooltip
    .html("<strong>Differences:</strong> <span style='color:red'>" + d.differences + "</span>")
    .style("left", (d3v3.mouse(this)[0]+70) + "px")
    .style("top", (d3v3.mouse(this)[1]) + "px")
}
var mouseleave = function(d) {
  Tooltip
    .style("opacity", 0)
  d3v3.select(this)
    .style("stroke", "none")
    .style("opacity", 0.8)
}



    var svg = d3v3.select("#degree_diff").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
    .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    //svg.call(tip);

    /*var fragment = document.createDocumentFragment();
    fragment.appendChild(document.getElementById('d3tip'));
    document.getElementById('degree_diff').appendChild(fragment);*/

    d3v3.tsv(path, type, function(error, data) {
    x.domain(data.map(function(d) { return d.nodes; }));
    y.domain([0, d3v3.max(data, function(d) { return d.differences; })]);

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis)
        .selectAll("text")
            .attr("transform", "rotate(90)")
            .attr("x", 10)
            .attr("y", -8)
            .style("text-anchor", "start");

        

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -40)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("Amount of differences");

    svg.selectAll(".bar")
        .data(data)
        .enter().append("rect")
        .attr("width", x.rangeBand())
        .attr("class", "bar")
        .attr("x", function(d) { return x(d.nodes); })
        .attr("y", function(d) { return y(d.differences); })
        .attr("height", function(d) { return height - y(d.differences); })
        .on("mouseover", mouseover)
        .on("mousemove", mousemove)
        .on("mouseleave", mouseleave)



    });

    function type(d) {
    d.differences = +d.differences;
    return d;
    }

}