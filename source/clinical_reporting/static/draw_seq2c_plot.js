
$(function() {
    var key = 'seq2c';

    function drawSeq2cPlot(data_el, placeholder_el, legend_placeholder_el) {
        var data = readJsonFromElement(data_el);

        //var nrmX = data.nrm.xs;
        //var nrmY = data.nrm.ys;
        //var ampX = data.amp.xs;
        //var ampY = data.amp.ys;
        //var ampG = data.amp.gs;
        //var delX = data.del.xs;
        //var delY = data.del.ys;
        //var delG = data.del.gs;

        var ticksX = data.ticksX;
        var linesX = data.linesX;
        var minY = data.minY;
        var maxY = data.maxY;

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

            for (var k = 0; k < data.events.length; k++) {
                var color = 'black';
                if (data.events[k].ampDel == 'Amp') color = 'green';
                if (data.events[k].ampDel == 'Del') color = 'red';

                var series = {
                    data: [[data.events[k].x, data.events[k].logRatio]],
                    geneName: data.events[k].geneName,
                    logRatio: data.events[k].logRatio,
                    ampDel: data.events[k].ampDel,
                    fragment: data.events[k].fragment,
                    color: color
                };
                series.points = {
                    show: true,
                    fill: true,
                    fillColor: color,
                    radius: 1
                };
                info.series.push(series);
            }

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
                        min: minY * 1.05,
                        max: maxY * 1.05,
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
                firstLabel.prepend('Log ratio' + '<span class="rhs">&nbsp;</span>=<span class="rhs">&nbsp;</span>');

                bindTip(placeholder_el, key, showTip, plot, 'top right', {});
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

    function showTip(item, plot, direction, generalData) {
        var LINE_HEIGHT = 16; // pixels

        direction = ((direction != null) ? direction : 'bottom right');
        //    pageY -= LINE_HEIGHT * (centralSeriesIndex + 1.5);
        var directions = direction.split(' ');
        // TODO: check screen width and tip width, if it doesn't fix - change position to "bottom left"

        var color = item.series.color;
        var geneName = item.series.geneName;
        var logRatio = item.series.logRatio;
        var ampDel = item.series.ampDel;
        var fragment = item.series.fragment;

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
                top: item.pageY + 5 - ((directions[0] == 'top') ? LINE_HEIGHT : 0),
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

        $('<div id="' + key + '_tip_line0">' +
            '<span style="z-index: 100; font-weight: bold; white-space: nowrap; color: ' + color + ';">' +
                geneName +
            '</span><br>' +
            '<span>' +
                'Log ratio = ' + toPrettyString(logRatio) +
                (ampDel ? ', ' + ampDel : '') +
                (fragment ? ', ' +  fragment : '') +
            '</span></div>')
            .css({
                height: LINE_HEIGHT})
            .appendTo('#' + key + '_plot_tip');
    }
    showTip.tipElementExists = false;

    drawSeq2cPlot($('#' + key + '_plot_data_json'),
                  $('#' + key + '_plot_placeholder'),
                  $('#' + key + '_plot_legend_placeholder'));
});

