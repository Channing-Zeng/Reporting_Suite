$(function() {
    var fullData = readJsonFromElement($('#substitutions_plot_data_json'));
    var colors = distinctColors();

    if (fullData != null) {
        if (fullData.hasOwnProperty('substitutions')) {
            drawSubstitutionsPlot('substitutions', fullData, $('#substitutions_plot_placeholder'));
        } else {
            var i = 0;
            for (var k in fullData) {
                if (fullData.hasOwnProperty(k)) {
                    var data = fullData[k];
                    drawSubstitutionsPlot(k + '_substitutions', data, $('#' + k + '_substitutions_plot_placeholder'));
                    $('#' + k + '_substitutions_header').css({color: colors[i]});
                    i += 1;
                }
            }
        }
    }
});

function drawSubstitutionsPlot(key, data, placeholder_el) {
    var minY = data.minY;
    var maxY = data.maxY;
    var lenTickX = 2.5;
    var barMaxWidth = 2;
    var subColors = ['#CC0000', '#CC6600', '#CCCC00', '#66CC00'];
    var maxRate = data.maxRate;

    var info = {
        isInitialized: false,
        maxY: 0,
        maxYTick: 0,
        series: null,
        showWithData: null
    };
    var substitutionsTypes = Object.keys(data.substitutions);
    var ticksX = [];
    substitutionsTypes.forEach(function(substitutionsType, i) {
        ticksX.push([lenTickX * i + barMaxWidth / 2, substitutionsType]);
    });
    var maxX = substitutionsTypes.length * lenTickX;
    if (!info.isInitialized) {
        info.series = [];

        substitutionsTypes.forEach(function(substitutionsType, k) {
            var color = subColors[Math.floor(k / 3)];

            var series = {
                data: [[lenTickX * k, data.substitutions[substitutionsType]]],
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
        });

        var series = {
          data: [[0, 0]],
          yaxis: 2
        };
        info.series.push(series);
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
                    mouseActiveRadius: 1000
                },
                yaxes: [{
                    position: "left",
                    min: minY,
                    max: maxY * 1.05,
                    labelWidth: 120,
                    reserveSpace: true,
                    lineWidth: 0.5,
                    color: '#000',
                    minTickSize: 1
                },
                {
                    position: "right",
                    tickFormatter: function (val, axis) {
                        return val + "%";
                    },
                    min: minY,
                    max: maxRate,
                    axisLabel: "Rate",
                    axisLabelUseCanvas: true,
                    axisLabelFontSizePixels: 12,
                    axisLabelPadding: 5
                }],
                xaxis: {
                    min: 0,
                    max: maxX,
                    tickLength: 1,
                    lineWidth: 0.5,
                    ticks: ticksX,
                    color: '#000'
                }
            });
            var firstLabel = placeholder_el.find('.y1Axis .tickLabel').last();
            firstLabel.prepend('Count' + '<span class="rhs">&nbsp;</span>=<span class="rhs">&nbsp;</span>');
        };
        info.isInitialized = true;
    }
    showPlotWithInfo(info);
}
