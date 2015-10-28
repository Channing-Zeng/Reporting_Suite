
$(function() {
    var key = 'gene';

    function drawGeneCovPlot(data_el, placeholder_el, legend_placeholder_el) {
        var data = readJsonFromElement(data_el);

        var geneNames = data.gene_names;

        var coordsX = data.coords_x;
        var aveDepths = data.ave_depths;
        var depthInThreshold = data.cov_in_thresh;
        var mutations = data.mutations;
        var cdsCovByGene = data.cds_cov_by_gene;

        var ticksX = data.ticks_x;
        var linesX = data.lines_x;

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

            for (var k = 0; k < coordsX.length; k++) {
                var percentDepthInThreshold = depthInThreshold[k] * 100;
                var curColor = getColorFromPercentCovered(percentDepthInThreshold);
                var series = {
                    data: [[coordsX[k], aveDepths[k]]],
                    label: geneNames[k],
                    color: curColor,
                    geneName: geneNames[k],
                    aveDepth: aveDepths[k],
                    mutations: mutations[geneNames[k]],
                    cdsDetails: cdsCovByGene[geneNames[k]],
                };
                series.points = {
                    show: true,
                    fill: true,
                    fillColor: curColor,
                    radius: 2
                };
                info.series.push(series);
            }

            var maxDepth = Math.max.apply(Math, aveDepths);
            info.showWithData = function(series, colors) {
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

                var firstLabel = placeholder_el.find('.yAxis .tickLabel').last();
                firstLabel.prepend('Ave depth' +
                    '<span class="rhs">&nbsp;</span>=<span class="rhs">&nbsp;</span>');

                bindTip(placeholder_el, key, showTip, plot, 'top right', {maxDepth: maxDepth});
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


    function getColorFromPercentCovered(percentCovered) {
        var hue = Math.max(0, (percentCovered - 50) * 2.4);
        var lightness = 33;
        var rgb = hslToRgb(hue / 360, 0.8, lightness / 100);
        return '#' + rgb[0].toString(16) + rgb[1].toString(16) + rgb[2].toString(16);
    }


    //data = {
    //    mutations: [],
    //    detailedCoverage: [{
    //        aveDepth: .0,
    //        percentInThreshold: .0
    //    }],
    //    maxDepth: 0
    //};


    function showTip(item, plot, direction, generalData) {
        var LINE_HEIGHT = 16; // pixels

        direction = ((direction != null) ? direction : 'bottom right');
        //    pageY -= LINE_HEIGHT * (centralSeriesIndex + 1.5);
        var directions = direction.split(' ');
        // TODO: check screen width and tip width, if it doesn't fix - change position to "bottom left"

        var maxDepth = generalData.maxDepth;
        var color = item.series.color;
        var geneName = item.series.geneName;
        var aveDepth = item.series.aveDepth;
        var cdsDetails = item.series.cdsDetails;
        var mutations = item.series.mutations;

        if (!showTip.tipElementExists) {
            $('<div id="' + key + '_plot_tip" class="white_stroked plot_tip"></div>')
                .appendTo('body')
                .css({'pointer-events': 'none'});

            $('<div id="' + key + '_plot_tip_vertical_rule" class="plot_tip_vertical_rule"></div>')
                .css({height: plot.height()})
                .appendTo('body');

            $('<div id="' + key + '_plot_tip_horizontal_rule" class="plot_tip_horizontal_rule"></div>')
                .css({width: plot.width()})
                .appendTo('body');

            showTip.tipElementExists = true;
        }

        $('#' + key + '_plot_tip')
            .html('')
            .css({
                top: item.pageY + 5 - ((directions[0] == 'top') ? LINE_HEIGHT * (mutations.length + 2) : 0),
                left: item.pageX + 10,
                zIndex: 1000})
            .show();

        $('#' + key + '_plot_tip_vertical_rule')
            .html('')
            .css({
                top: plot.offset().top,
                left: item.pageX})
            .show();

        $('#' + key + '_plot_tip_horizontal_rule')
            .html('')
            .css({
                top: item.pageY,
                left: plot.offset().left})
            .show();

        $('<div id="' + key + '_tip_line0"><span style="z-index: 100; white-space: nowrap; ' +
            'color: ' + item.color + ';">' + geneName + '</span> - ' +
            toPrettyString(aveDepth, 'x') + ' ave depth</div>')
            .css({
                height: LINE_HEIGHT,
                'font-weight': 'bold'})
            .appendTo('#' + key + '_plot_tip');

        for (var m = 0; m < mutations.length; m++) {
            var mut = mutations[m];
            if (mut != '.') {
                $('<div id="' + key + '_tip_line' + m + '"><span>' + mut + '</span></div>')
                    .css({
                        height: LINE_HEIGHT,
                        'white-space': 'nowrap'})
                    .appendTo('#' + key + '_plot_tip');
            }
        }
        //(plot show' + '0' +
        //    '...' + toPrettyString(maxDepth) + 'x)
        $('<div id="' + key + '_tip_line' + (m + 1) + '">' +
            '<span>CDS coverage across the gene ' +
            '(total ' + cdsDetails.length + ' CDS)</span></div>')
            .css({
                height: LINE_HEIGHT,
                'white-space': 'nowrap',
                'margin-top': '5px'})
            .appendTo('#' + key + '_plot_tip');

        var miniSeries = [];
        var prevX = 0;

        for (var k = 0; k < cdsDetails.length; k++) {
            var percentInThreshold = cdsDetails[k].percentInThreshold * 100;
            var curColor = getColorFromPercentCovered(percentInThreshold);
            var cdsAveDepth = cdsDetails[k].aveDepth;
            var cdsSize = cdsDetails[k].end - cdsDetails[k].start;
            //minY = Math.min(y, minY);
            //maxY = Math.max(y, maxY);
            var bar = {
                data: [
                    [prevX, cdsAveDepth],
                    [prevX + cdsSize, cdsAveDepth]
                ],
                color: curColor
            };
            bar.lines = {
                show: true,
                lineWidth: 1
            };
            prevX = prevX + cdsSize;
            miniSeries.push(bar);
        }
        //console.log('Total CDS: ' + cdsDetails.length + ', total size: ' + prevX);

        var miniPlotWidth = prevX / 15;
        $('<div class="mini_plot_placeholder" id="mini_' + key + '_plot_placeholder"></div>')
            .appendTo('#' + key + '_plot_tip')
            .css({
                width: miniPlotWidth + 'px',
                height: '60px',
                'margin-left': '-5px'});

        $.plot($('#mini_' + key + '_plot_placeholder'), miniSeries, {
            shadowSize: 0,
            grid: {
                borderWidth: 1,
                hoverable: true,
                autoHighlight: false,
                mouseActiveRadius: 1000,
                backgroundColor: '#FFF',
            },
            yaxis: {
                min: 0,
                max: maxDepth + 1,
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

    showTip.tipElementExists = false;

    drawGeneCovPlot($('#' + key + '_plot_data_json'),
                    $('#' + key + '_plot_placeholder'),
                    $('#' + key + '_plot_legend_placeholder'))
});
