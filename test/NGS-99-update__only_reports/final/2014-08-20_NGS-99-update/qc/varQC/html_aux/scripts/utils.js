
/**********/
/* COLORS */

// var colors = ["#FF5900", "#008FFF", "#168A16", "#7C00FF", "#00B7FF", "#FF0080", "#7AE01B", "#782400", "#E01B6A"];
var colors = [
    '#FF0000', //red
    '#0000FF', //blue
    '#008000', //green
    '#FFA500', //orange
    '#FF00FF', //fushua
    '#CCCC00', //yellow
    '#800000', //maroon
    '#00CCCC', //aqua
    '#808080', //gray
    '#800080', //purple
    '#808000', //olive
    '#000080', //navy
    '#008080', //team
    '#00FF00', //lime
];

function distinctColors(count) {
    var colors = [];
    for(var hue = 0; hue < 360; hue += 360 / count) {
        var color = hsvToRgb(hue, 100, 100);
        var colorStr = '#' + color[0].toString(16) + color[1].toString(16) + color[2].toString(16);
        colors.push();
    }
    return colors;
}

/**************/
/* FORMATTING */
function isIntegral(num) {
    return num % 1 === 0;
}

function isFractional(num) {
    return !isIntegral(num);
}

function toPrettyString(num, unit) {
    var str,
        frac_digits = 0;

    if (typeof num === 'number') {
        if (num <= 999) {
            if (isFractional(num)) {
                if (isIntegral(num * 10)) {
                    frac_digits = 1;

                } else if (isIntegral(num * 100) || num > 100) {
                    frac_digits = 2;

                } else {
                    var lessthen = 0.00000000001;
                    frac_digits = 12;

                    while (num > lessthen && frac_digits > 0) {
                        frac_digits--;
                        lessthen *= 10;
                    }
                    if (!isIntegral(num * Math.pow(10, frac_digits))) {
                        frac_digits++;
                        if (!isIntegral(num * Math.pow(10, frac_digits))) {
                            frac_digits++;
                        }
                    }
                }
            }
            str = num.toFixed(frac_digits);
        } else {
            str = num.toFixed(0).replace(/(\d)(?=(\d\d\d)+(?!\d))/g,'$1<span class=\'hs\'></span>');
        }
        str += (unit ? '<span class=\'rhs\'>&nbsp;</span>' + unit : '');
    } else {
        str = num;
    }
    return str;
}

function ordinalNumberToPrettyString(num, unit) {
    var numStr = num.toString();
    var lastDigit = numStr[numStr.length-1];
    var beforeLastDigit = numStr[numStr.length-2];

    var res = toPrettyString(num);

    if (lastDigit == '1' && beforeLastDigit != '1') {
        res += "st";
    } else if (lastDigit == '2' && beforeLastDigit != '1') {
        res += "nd";
    } else if (lastDigit == '3' && beforeLastDigit != '1') {
        res += "rd";
    } else {
        res += 'th';
    }

    res += (unit ? '<span class="rhs">&nbsp;</span>' + unit : '');

    return res;
}

function getMaxDecimalTick(maxY) {
    var maxYTick = maxY;
    if (maxY <= 100000000000) {
        maxYTick = Math.ceil((maxY+1)/10000000000)*10000000000;
    } if (maxY <= 10000000000) {
        maxYTick = Math.ceil((maxY+1)/1000000000)*1000000000;
    } if (maxY <= 1000000000) {
        maxYTick = Math.ceil((maxY+1)/100000000)*100000000;
    } if (maxY <= 100000000) {
        maxYTick = Math.ceil((maxY+1)/10000000)*10000000;
    } if (maxY <= 10000000) {
        maxYTick = Math.ceil((maxY+1)/1000000)*1000000;
    } if (maxY <= 1000000) {
        maxYTick = Math.ceil((maxY+1)/100000)*100000;
    } if (maxY <= 100000) {
        maxYTick = Math.ceil((maxY+1)/10000)*10000;
    } if (maxY <= 10000) {
        maxYTick = Math.ceil((maxY+1)/1000)*1000;
    } if (maxY <= 1000) {
        maxYTick = Math.ceil((maxY+1)/100)*100.
    } if (maxY <= 100) {
        maxYTick = Math.ceil((maxY+1)/10)*10.
    }
    return maxYTick;
}

function getBpTickFormatter(maxY, additionalText) {
    additionalText = additionalText || '';

    return function(val, axis) {
        var res;

        if (val == 0) {
            res = 0;

        } else if (val >= 1000000) {
            res = val / 1000000;

            if (val > maxY + 1 || val + axis.tickSize >= 1000000000) {
                res = additionalText + toPrettyString(res, 'Mbp');
            } else {
                res = toPrettyString(res);
            }
        } else if (val >= 1000) {
            res = val / 1000;

            if (val > maxY + 1 || val + axis.tickSize >= 1000000) {
                res = additionalText + toPrettyString(res, 'kbp');
            } else {
                res = toPrettyString(res);
            }
        } else if (val >= 1) {
            res = val;

            if (val > maxY + 1 || val + axis.tickSize >= 1000) {
                res = additionalText + toPrettyString(res, 'bp');
            } else {
                res = toPrettyString(res);
            }
        }
        return res;
    }
}

function windowsTickFormatter(v, axis) {
    return toPrettyString(v);
//    var val = v.toFixed(0);
//    if (!gc.yAxisLabeled && val > gc.maxY) {
//        gc.yAxisLabeled = true;
//        var res = val + ' window';
//        if (val > 1) {
//            res += 's'
//        }
//        return res;
//    } else {
//        return val;
//    }
}

function getBpLogTickFormatter(maxY) {
    return getBpTickFormatter(maxY);
}

function getContigNumberTickFormatter(maxX) {
    return function (val, axis) {
        if (typeof axis.tickSize == 'number' && val > maxX - axis.tickSize) {
            return "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + ordinalNumberToPrettyString(val, 'contig');
        } else {
            return val;
        }
    }
}

function trim(str) {
    return str.replace(/^\s+/g, '');
}

function nbsp(str, metricName) {
    if (metricName.length > 0 && metricName[0] == ' ') {
        str = '&nbsp;&nbsp;&nbsp;' + str;
    }
    return str;
}

function containsObject(obj, list) {
    var i;
    for (i = 0; i < list.length; i++) {
        if (list[i] === obj) {
            return true;
        }
    }

    return false;
}

/*********************/
/* GLOSSARY TOOLTIPS */
function addTooltipIfDefinitionExists(glossary, metricName) {
    metricName = trim(metricName);

    if (containsObject(metricName, Object.keys(glossary))) {
        return '<a class="tooltip-link" rel="tooltip" title="' +
            metricName + ' ' + glossary[metricName] + '">' + metricName + '</a>';
    } else {
        return metricName;
    }
}

/*************/
/* PLOT TIPS */
function bindTip(placeholder, series, plot, xToPrettyStringFunction, xUnit, position) {
    var prevPoint = null;

    $(placeholder).bind("plothover", function(event, pos, item) {
        if (dragTable && dragTable.isDragging)
            return;

        if (item) {
            if (prevPoint != item.dataIndex) {
                prevPoint = item.dataIndex;

                var x = item.datapoint[0];

                showTip(item.pageX, item.pageY, plot.offset(),
                    plot.width(), plot.height(),
                    series, item.seriesIndex, item.dataIndex,
                    xToPrettyStringFunction(x, xUnit) + ':',
                    position);
            }
        } else {
            $('#plot_tip').hide();
            $('#plot_tip_vertical_rule').hide();
            $('#plot_tip_horizontal_rule').hide();
            prevPoint = null;
        }
    });
}

var tipElementExists = false;
function showTip(pageX, pageY, offset, plotWidth, plotHeight,
                 series, centralSeriesIndex, xIndex, xStr, position) {
    var LINE_HEIGHT = 16; // pixels

    position = ((position != null) ? position : 'bottom right');
//    pageY -= LINE_HEIGHT * (centralSeriesIndex + 1.5);

    var directions = position.split(' ');

    if (!tipElementExists) {
        $('<div id="plot_tip" class="white_stroked"></div>').appendTo('body');

        $('<div id="plot_tip_vertical_rule"></div>').css({
            height: plotHeight
        }).appendTo('body');

        $('<div id="plot_tip_horizontal_rule"></div>').css({
            width: plotWidth
        }).appendTo('body');

        tipElementExists = true;
    }

    $('#plot_tip').html('').css({
        top: pageY + 5 - ((directions[0] == 'top') ? LINE_HEIGHT * (series.length + 2) : 0),
        left: pageX + 10
    }).show();

    $('#plot_tip_vertical_rule').html('').css({
        top: offset.top,
        left: pageX
    }).show();

    $('#plot_tip_horizontal_rule').html('').css({
        top: pageY,
        left: offset.left
    }).show();

    $('<div>' + xStr + '</div>').css({
        height: LINE_HEIGHT
    }).appendTo('#plot_tip');

    var sortedYsAndColors = [];
    for (var i = 0; i < series.length; i++) {
        sortedYsAndColors.push({
            y: (series[i].data[xIndex] || series[i].data[series[i].data.length - 1])[1],
            color: series[i].color,
            label: (series[i].isReference ? 'Reference' : series[i].label),
            isCurrent: i == centralSeriesIndex
        });
    }
    sortedYsAndColors.sort(function(a, b) { return a.y < b.y;});

    for (i = 0; i < sortedYsAndColors.length; i++) {
        var item =sortedYsAndColors[i];

        $('<div id="tip_line' + i + '">' + toPrettyString(item.y)
            + ', <span style="color: ' + item.color + ';">' + item.label + '</span></div>').css({
                height: LINE_HEIGHT,
                "font-weight": item.isCurrent ? "bold" : "normal"
            }).appendTo('#plot_tip');
    }
}



// Cookie functions based on http://www.quirksmode.org/js/cookies.html
// Cookies won't work for local files.

var createCookie = function(name, value, days) {
    var expires = '';
    if (days) {
        var date = new Date();
        date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
        expires = '; expires=' + date.toGMTString();
    }
    var path = document.location.pathname;
    document.cookie = name + '=' + value + expires + '; path=' + path;
};

var readCookie = function(name) {
    var nameEQ = name + '=';
    var ca = document.cookie.split(';');
    for(var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) == ' ') {
            c = c.substring(1, c.length);
        }
        if (c.indexOf(nameEQ) == 0) {
            return c.substring(nameEQ.length, c.length);
        }
    }
    return null;
};

var eraseCookie = function(name) {
    createCookie(name, '', -1);
};


// Dean's forEach: http://dean.edwards.name/base/forEach.js
/*forEach, version 1.0
 Copyright 2006, Dean Edwards
 License: http://www.opensource.org/licenses/mit-license.php */

// array-like enumeration
if (!Array.forEach) { // mozilla already supports this
    Array.forEach = function(array, block, context) {
        for (var i = 0; i < array.length; i++) {
            block.call(context, array[i], i, array);
        }
    };
}

// generic enumeration
Function.prototype.forEach = function(object, block, context) {
    for (var key in object) {
        if (typeof this.prototype[key] == "undefined") {
            block.call(context, object[key], key, object);
        }
    }
};

// character enumeration
String.forEach = function(string, block, context) {
    Array.forEach(string.split(""), function(chr, index) {
        block.call(context, chr, index, string);
    });
};

// globally resolve forEach enumeration
var forEach = function(object, block, context) {
    if (object) {
        if (object instanceof Function) {                 // functions have a "length" property
            Function.forEach(object, block, context);

        } else if (object.each instanceof Function) {     // jQuery
            object.each(function(i, elt) {
                block(elt);
            }, context);

        } else if (object.forEach instanceof Function) {  // the object implements a custom forEach method
            object.forEach(block, context);

        } else if (typeof object == "string") {           // a string
            String.forEach(object, block, context);

        } else if (typeof object.length == "number") {    // array-like object
            Array.forEach(object, block, context);

        } else {
            Object.forEach(object, block, context);
        }
    }
};


jQuery.fn.exists = function(){
    return jQuery(this).length > 0;
};


// http://krijnhoetmer.nl/stuff/javascript/table-align-char/
//function make_table_align_char() {
//     var currencies = /(\$|€|&euro;)/;
//     var leftWidth = 0, rightWidth = 0, currencyWidth;
//     for (var tableCounter = 0, tables = document.getElementsByTagName('table'); tableCounter < tables.length; tableCounter++) {
//      if (tables[tableCounter].className.indexOf('fix-align-char') != -1) {
//       var fCols = [];
//       for (var i = 0, cols = tables[tableCounter].getElementsByTagName('col'); i < cols.length; i++) {
//        if (cols[i].getAttribute('char')) {
//         fCols[i] = cols[i].getAttribute('char');
//        };
//       };
//       var leftPart, rightPart, parts;
//       for (var i = 0, trs = tables[tableCounter].rows; i < trs.length; i++) {
//        for (var j = 0, tds = trs[i].getElementsByTagName('td'); j < tds.length; j++) {
//         if (fCols[j]) {
//          if (tds[j].innerHTML.indexOf(fCols[j]) != -1) {
//           parts = tds[j].innerHTML.split(fCols[j]);
//           leftPart = parts.slice(0, parts.length -1).join(fCols[j]);
//           leftPart = leftPart.replace(currencies, '<span class="currency">$1</span>');
//           rightPart = fCols[j] + parts.pop();
//           tds[j].innerHTML = '<span class="left">' + leftPart + '</span><span class="right">' + rightPart + '</span>';
//          } else {
//           tds[j].innerHTML = tds[j].innerHTML.replace(currencies, '<span class="currency">$1</span>');
//           tds[j].innerHTML = '<span class="left">' + tds[j].innerHTML + '</span>';
//          };
//          tds[j].className = 'char-align';
//          var txt = document.createTextNode(tds[j].firstChild.offsetWidth);
//          if (leftWidth < tds[j].firstChild.offsetWidth) {
//           leftWidth = tds[j].firstChild.offsetWidth;
//          };
//          if (tds[j].childNodes[1]) {
//           var txt = document.createTextNode(tds[j].childNodes[1].offsetWidth);
//           if (rightWidth < tds[j].childNodes[1].offsetWidth) {
//            rightWidth = tds[j].childNodes[1].offsetWidth;
//           };
//          };
//         };
//        };
//       };
//      };
//     };
//     // This is ugly and should be improved (amongst other parts of the code ;)
//     var styleText = '<style type="text/css">.fix-align-char td.char-align { width: ' + (leftWidth + rightWidth) + 'px; }\n.fix-align-char span.left { float: left; text-align: right; width: ' + (leftWidth) + 'px; }\n.fix-align-char span.currency { text-align: left; float: left; }\n.fix-align-char span.right { float: right; text-align: left; width: ' + rightWidth + 'px; }</style>';
//     document.body.innerHTML += styleText;
//}
//
//
//


