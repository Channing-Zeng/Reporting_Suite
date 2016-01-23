$(function() {
    var fullData = readJsonFromElement($('#mut_plot_data_json'));
    var colors = distinctColors();

    if (fullData != null) {
        if (fullData.hasOwnProperty('mutations')) {
            drawMutPlot('mut', fullData, $('#mut_plot_placeholder'));
        } else {
            var i = 0;
            for (var k in fullData) {
                if (fullData.hasOwnProperty(k)) {
                    var data = fullData[k];
                    drawMutPlot(k + '_mut', data, $('#' + k + '_mut_plot_placeholder'));
                    $('#' + k + '_mut_header').css({color: colors[i]});
                    i += 1;
                }
            }
        }
    }
});

function drawMutPlot(key, data, placeholder_el) {
    var minY = data.minY;
    var maxY = 100;
    var maxX = 100;
    var ticksX = maxX / data.mutations.length;
    var barMaxWidth = ticksX / 1.15;

    var info = {
        isInitialized: false,
        maxY: 0,
        maxYTick: 0,
        series: null,
        showWithData: null
    };
    var mutationTypes = ['frameshift variant', 'inframe variant', 'nonsense variant', 'missense variant', 'stop gained',
        'splice site', 'start lost', 'UTR variant'];
    mutationTypes.forEach(function(mutType, i) {
        var id = 'label_' + i + '_id';
        $('#mut_plot_legend_placeholder').append('<div style="display: inline-block; margin-right: 2em;">' +
            '<div style="display: inline-block; width: 5px; height: 5px; border-radius: 5px; background-color: ' + colors[i] + '"></div>' +
            '&nbsp<label for="' + id + '" style="color: ' + colors[i] + '">' + mutType  + '</label>' +
            '</div>');
    });
    if (!info.isInitialized) {
        info.series = [];

        for (var k = 0; k < data.mutations.length; k++) {
            var color = 'black';

            mutationTypes.forEach(function(mutType, i) {
                if (data.mutations[k].mutType.indexOf(mutType.split(' ')[0]) != -1)
                color = colors[i];
            });

            var series = {
                data: [[ticksX * k, data.mutations[k].freq]],
                geneName: data.mutations[k].geneName,
                freq: data.mutations[k].freq,
                position: data.mutations[k].position,
                mutType: data.mutations[k].mutType,
                aaChg: data.mutations[k].aaChg,
                cdnaChange: data.mutations[k].cdnaChange,
                color: color
            };
            series.bars = {
                show: true,
                barWidth: barMaxWidth,
                lineWidth: 0,
                order: 1,
                fillColor: color
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
                    autoHighlight: false
                },
                yaxis: {
                    min: minY,
                    max: maxY,
                    labelWidth: 120,
                    reserveSpace: true,
                    lineWidth: 0.5,
                    color: '#000',
                    minTickSize: 1
                },
                xaxis: {
                    min: 0,
                    max: maxX,
                    tickLength: 1,
                    lineWidth: 0.5,
                    color: '#000'
                }
            });
            var firstLabel = placeholder_el.find('.yAxis .tickLabel').last();
            firstLabel.prepend('Frequency' + '<span class="rhs">&nbsp;</span>=<span class="rhs">&nbsp;</span>');
            firstLabel.append('%');

            bindTip(placeholder_el, key, showMutationTip, plot, 'top right', {key: key});
        };

        info.isInitialized = true;
    }

    showPlotWithInfo(info);

    var subt_placeholder = $('#substitutions_plot_placeholder');
    var subst_plot_width = subt_placeholder.width();
    var mut_plot_width = placeholder_el.width();
    var mut_plot_pos = placeholder_el.offset();
    console.log(subst_plot_width);
    console.log(mut_plot_pos);
    mut_plot_pos.left += mut_plot_width - subst_plot_width + 151;
    subt_placeholder
        .css('position', 'absolute')
        .css(mut_plot_pos);
        //.css('right', mut_plot_right - subst_plot_width);
}

function showMutationTip(key, item, plot, direction, generalData) {
    var LINE_HEIGHT = 16; // pixels

    direction = ((direction != null) ? direction : 'bottom right');
    //    pageY -= LINE_HEIGHT * (centralSeriesIndex + 1.5);
    var directions = direction.split(' ');
    // TODO: check screen width and tip width, if it doesn't fix - change position to "bottom left"

    var color = item.series.color;
    var geneName = item.series.geneName;
    var frequency = item.series.freq;
    var change = item.series.aaChg;
    var cdnaChange = item.series.cdnaChange;

    if (!showMutationTip.tipElementExists) {
    $('<div id="' + key + '_plot_tip" class="white_stroked plot_tip"></div>')
        .appendTo('body')
        .css({'pointer-events': 'none'});

    $('<div id="' + key + '_plot_tip_vertical_rule" class="plot_tip_vertical_rule"></div>')
        .css({height: plot.height()})
        .appendTo('body');

    $('<div id="' + key + '_plot_tip_horizontal_rule" class="plot_tip_horizontal_rule"></div>')
        .css({width: plot.width()})
        .appendTo('body');
    showMutationTip.tipElementExists = true;
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
            geneName + '</span> - ' + toPrettyString(frequency) + '% AF</div>')
        .css({
            height: LINE_HEIGHT})
        .appendTo('#' + key + '_plot_tip');
    $('<div id="' + key + '_tip_line1">' +
        '<span>' + (change ? 'p.' + change : cdnaChange) + '</span></div>')
        .css({
            height: LINE_HEIGHT})
        .appendTo('#' + key + '_plot_tip');
}

