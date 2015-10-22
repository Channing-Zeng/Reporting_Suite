$(function() {
    drawGeneCovPlot($('#gene_plot_data_json'), $('#gene_plot_placeholder'), $('#gene_plot_legend_placeholder'))
});

function drawGeneCovPlot(data_el, placeholder_el, legend_placeholder_el) {
    var data = readJsonFromElement(data_el);

    var geneNames = data.gene_names;

    var coordX = data.coord_x;
    var coordY = data.coord_y;
    var depthInThreshold = data.cov_in_thresh;

    var ticksX = data.ticks_x;
    var linesX = data.lines_x;

    var mutations = data.mutations;
    var detailedCoverage = data.detailed_cov;
    var info = {
        isInitialized: false,
        maxY: 0,
        maxYTick: 0,
        series: null,
        showWithData: null
    };

    var gridLines = [];
    for (var i = 0; i < linesX.length; i++) {
        var line = { xaxis: { from: linesX[i], to: linesX[i] }, color: "#ccc", lineWidth: 0.5 };
        gridLines.push(line);
    }

    if (!info.isInitialized) {
        info.series = [];

        for (var k = 0; k < coordX.length; k++) {
            var percentDepthInThreshold = depthInThreshold[k] * 100;
            var curColor = getColorFromPercentCovered(percentDepthInThreshold);
            var series = {
                data: [[coordX[k], coordY[k]]],
                label: geneNames[k],
                color: curColor
            };
            series.points = {
                show: true,
                fill: true,
                fillColor: curColor,
                radius: 2
            };
            info.series.push(series);
        }

        var maxDepth = Math.max.apply(Math, coordY);
        info.showWithData = function (series, colors) {
            var plot = $.plot(placeholder_el, series, {
                shadowSize: 0,
                colors: colors,
                legend: {
                    container: $('useless-invisible-element-that-does-not-even-exist')
                },
                grid: {
                    borderWidth: 1,
                    hoverable: true,
                    autoHighlight: false,
                    mouseActiveRadius: 1000,
                    markings: gridLines
                },
                yaxis: {
                    min: 0,
                    max: maxDepth * 1.1,
                    labelWidth: 120,
                    reserveSpace: true,
                    lineWidth: 0.5,
                    color: '#000',
                    minTickSize: 1
                },
                xaxis: {
                    min: 0,
                    max: Math.max.apply(Math, linesX),
                    ticks: ticksX,
                    tickLength: 0,
                    lineWidth: 0.5,
                    color: '#000'
                }
            });

            var firstLabel = $('.yAxis .tickLabel').last();
            firstLabel.prepend('Ave depth' + '<span class="rhs">&nbsp;</span>=<span class="rhs">&nbsp;</span>');

            bindTip(placeholder_el, series, mutations, detailedCoverage, plot, toPrettyString, 1, '', 'top right', maxDepth);
        };

        info.isInitialized = true;
    }

    $.each(info.series, function(i, series) {
        legend_placeholder_el.find('#label_' + series.number + '_id').click(function() {
            showPlotWithInfo(info);
        });
    });

    showPlotWithInfo(info);
}

function showPlotWithInfo(info, index) {
    var series = info.series;
    var colors = [];
    for (var i = 0; i < series.length; i++) {
        colors.push(series[i].color);
    }
    info.showWithData(info.series, colors);

}

function getColorFromPercentCovered(percentCovered) {
    var hue = Math.max(0, (percentCovered - 50) * 2.4);
    var lightness = 33;
    var rgb = hslToRgb(hue / 360, 0.8, lightness / 100);
    return '#' + rgb[0].toString(16) + rgb[1].toString(16) + rgb[2].toString(16);
}

function hslToRgb(h, s, l){
    var r, g, b;

    if(s == 0) {
        r = g = b = l; // achromatic
    } else {
        function hue2rgb(p, q, t){
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        }

        var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        var p = 2 * l - q;
        r = hue2rgb(p, q, h + 1/3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1/3);
    }

    return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
}

function readJsonFromElement(element) {
    var result;
    try {
        result = JSON.parse($(element).html());
    } catch (e) {
        result = null;
    }
    return result;
}

function bindTip(placeholder, series, mutations, detailedCoverage, plot, xToPrettyStringFunction, tickX, xUnit, position, maxDepth) {
    var prevPoint = null;

    $(placeholder).bind("plothover", function(event, pos, item) {
        if (item) {
          if (prevPoint != item.seriesIndex) {
              prevPoint = item.seriesIndex;
              var x = item.datapoint[0];
              showTip(item.pageX, item.pageY, plot.offset(),
                  plot.width(), plot.height(),
                  series, item.seriesIndex, x, item.dataIndex,
                  xToPrettyStringFunction(x, xUnit, tickX) + ':',
                  position, mutations, detailedCoverage, maxDepth);
            }
        } else {
            $('#plot_tip').hide();
            $('#plot_tip_vertical_rule').hide();
            $('#plot_tip_horizontal_rule').hide();
            prevPoint = null;
        }
    });
}

function showTip(pageX, pageY, offset, plotWidth, plotHeight,
    series, centralSeriesIndex, xPos, xIndex, xStr, position, mutations, detailedCoverage, maxDepth) {
    var LINE_HEIGHT = 16; // pixels

    position = ((position != null) ? position : 'bottom right');
    //    pageY -= LINE_HEIGHT * (centralSeriesIndex + 1.5);

    var directions = position.split(' ');

    var item = {
        aveDepth: series[centralSeriesIndex].data[xIndex][1],
        color: series[centralSeriesIndex].color,
        label: series[centralSeriesIndex].label
    };
    var geneName = item.label;
    var geneMutations = mutations[geneName];
    if (!geneMutations) geneMutations = [];

    if (!tipElementExists) {
        $('<div id="plot_tip" class="white_stroked"></div>').appendTo('body');

        $('<div id="plot_tip_vertical_rule"></div>').css({
            height: plotHeight,
        }).appendTo('body');

        $('<div id="plot_tip_horizontal_rule"></div>').css({
            width: plotWidth,
        }).appendTo('body');

        tipElementExists = true;
    }

    $('#plot_tip').html('').css({
        top: pageY + 5 - ((directions[0] == 'top') ? LINE_HEIGHT * (geneMutations.length + 2) : 0),
        left: pageX + 10,
        zIndex: 1000
    }).show();

    $('#plot_tip_vertical_rule').html('').css({
        top: offset.top,
        left: pageX,
    }).show();

    $('#plot_tip_horizontal_rule').html('').css({
        top: pageY,
        left: offset.left,
    }).show();
    /*
    $('<div>' + xStr + '</div>').css({
        height: LINE_HEIGHT,
    }).appendTo('#plot_tip');*/

    var newEl = '<div id="tip_line0"><span style="z-index: 100; color: ' + item.color + ';">' + geneName + '</span> - ' + toPrettyString(item.aveDepth) + 'x ave depth</div>';
    console.log(newEl);
    $(newEl).css({
        height: LINE_HEIGHT,
        "font-weight": "bold"
    }).appendTo('#plot_tip');
    for (var m = 0; m < geneMutations.length; m++) {
        $('<div id="tip_line"' + m + '><span>' + geneMutations[m].join(', ') + '</span></div>').css({
          height: LINE_HEIGHT
        }).appendTo('#plot_tip');
    }
    $('<div id="tip_line"' + (m + 1) + ' style="margin-top: 5px;"><span>CDS coverage (' + '0' + '...' + toPrettyString(maxDepth) + 'x) accross the gene</span></div>').css({
        height: LINE_HEIGHT
    }).appendTo('#plot_tip');
    var geneCoverage = detailedCoverage[geneName];
    var miniSeries = [];
    var minY = 0;
    var maxY = maxDepth;
    console.log(maxDepth);
    var prevX = 0;

    for (var k = 0; k < geneCoverage.length; k++) {
        var percentageInThreshold = geneCoverage[k][2] * 100;
        var curColor = getColorFromPercentCovered(percentageInThreshold);
        var aveDepth = geneCoverage[k][1];
        //minY = Math.min(y, minY);
        //maxY = Math.max(y, maxY);
        var bar = {
            data: [[prevX, geneCoverage[k][1]], [prevX + geneCoverage[k][0], aveDepth]],
            color: curColor
        };
        bar.lines = {
            show: true,
            lineWidth: 1
        };
        prevX = prevX + geneCoverage[k][0];
        miniSeries.push(bar);
    }

    var geneLength = prevX;
    var miniPlotWidth = Math.max(geneLength / 15);
    $('<div class="mini_plot_placeholder" ' +
        'id="mini_gene_plot_placeholder" ' +
        'style="width:' + miniPlotWidth + 'px; height: 60px; margin-left: -5px;"></div>').appendTo('#plot_tip');

    $.plot($("#mini_gene_plot_placeholder"), miniSeries, {
        shadowSize: 0,
        grid: {
            borderWidth: 1,
            hoverable: true,
            autoHighlight: false,
            mouseActiveRadius: 1000,
            backgroundColor: '#FFF',
        },
        yaxis: {
            min: minY * 0.9,
            max: maxY * 1.1,
            labelWidth: 0,
            reserveSpace: true,
            ticks: [],
            lineWidth: 0.5,
            minTickSize: 1
        },
        xaxis: {
            min: 0,
            max: prevX,
            ticks: [],
            tickLength: 0,
            lineWidth: 0.5
        }
    });
}

