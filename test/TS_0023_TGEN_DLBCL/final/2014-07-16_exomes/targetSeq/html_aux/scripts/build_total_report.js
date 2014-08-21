// Generated by CoffeeScript 1.7.1
(function() {
  var CSS_PROP_TO_COLOR, DRAGGABLE_COLUMNS, GREEN_HSL, GREEN_HUE, RED_HUE, calc_cell_contents, calc_records_cell_contents, get_color, get_meta_tag_contents, get_metric_name_html, metric, record, report, set_legend,
    __hasProp = {}.hasOwnProperty;

  report = {
    sample: {
      name: '',
      phenotype: '',
      bam: '',
      bed: '',
      vcf_by_caller: {
        name: '',
        summary_qc_rep_fpaths: [],
        anno_vcf_fpaths: {},
        anno_filt_vcf_fpaths: {}
      }
    },
    fpath: '',
    link: '',
    records: []
  };

  record = {
    metric: null,
    value: '',
    meta: null
  };

  metric = {
    name: '',
    short_name: '',
    description: '',
    quality: '',
    presision: 0,
    type: null
  };

  DRAGGABLE_COLUMNS = true;

  RED_HUE = 0;

  GREEN_HUE = 120;

  GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)';

  CSS_PROP_TO_COLOR = 'background-color';

  get_color = function(hue) {
    var lightness;
    lightness = 92;
    return 'hsl(' + hue + ', 80%, ' + lightness + '%)';
  };

  get_meta_tag_contents = function(rec) {
    var a, db, dbs, k, meta, meta_table, novelty, val, values, _i, _len;
    meta = rec.meta;
    if ((meta != null) && ((function() {
      var _results;
      _results = [];
      for (a in meta) {
        _results.push(a);
      }
      return _results;
    })()).length !== 0) {
      if (typeof meta === 'string') {
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"" + meta + "\"";
      } else {
        ((function() {
          var _results;
          _results = [];
          for (k in meta) {
            if (!__hasProp.call(meta, k)) continue;
            _results.push(k);
          }
          return _results;
        })()).length !== 0;
        meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>';
        for (novelty in meta) {
          values = meta[novelty];
          if (novelty !== 'all') {
            meta_table += "<td>" + novelty + "</td>";
          }
        }
        meta_table += '</tr>\n';
        for (novelty in meta) {
          values = meta[novelty];
          dbs = (function() {
            var _results;
            _results = [];
            for (db in values) {
              val = values[db];
              if (db !== 'average') {
                _results.push(db);
              }
            }
            return _results;
          })();
          dbs.push('average');
          break;
        }
        for (_i = 0, _len = dbs.length; _i < _len; _i++) {
          db = dbs[_i];
          meta_table += "<tr><td>" + db + "</td>";
          for (novelty in meta) {
            values = meta[novelty];
            if (novelty !== 'all') {
              meta_table += "<td>" + (toPrettyString(values[db], rec.metric.unit)) + "</td>";
            }
          }
          meta_table += '</tr>\n';
        }
        meta_table += '</table>\n';
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"" + meta_table + "\"";
      }
    } else {
      return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\"";
    }
  };

  get_metric_name_html = function(rec, use_full_name) {
    var description, metricName;
    if (use_full_name == null) {
      use_full_name = false;
    }
    if (rec.metric.short_name && !use_full_name) {
      metricName = rec.metric.short_name;
      description = rec.metric.description || rec.metric.name;
      return "<a class=\"tooltip-link\" rel=\"tooltip\" title=\"" + description + "\">" + metricName + "</a>";
    } else {
      return rec.metric.name;
    }
  };

  calc_records_cell_contents = function(records, font) {
    var num_html, rec, result, value, _i, _len, _results;
    _results = [];
    for (_i = 0, _len = records.length; _i < _len; _i++) {
      rec = records[_i];
      value = rec.value;
      num_html = '';
      if ((value == null) || value === '') {
        rec.cell_contents = '-';
      } else {
        if (typeof value === 'number') {
          rec.num = value;
          rec.cell_contents = toPrettyString(value, rec.metric.unit);
          num_html = toPrettyString(value);
        } else if (/^-?.?[0-9]/.test(value)) {
          result = /([0-9\.]+)(.*)/.exec(value);
          rec.num = parseFloat(result[1]);
          rec.cell_contents = toPrettyString(rec.num, rec.metric.unit) + result[2];
          num_html = toPrettyString(rec.num);
        } else {
          rec.cell_contents = value;
        }
      }
      _results.push(rec.frac_width = $.fn.intPartTextWidth(num_html, font));
    }
    return _results;
  };

  calc_cell_contents = function(report, font) {
    var a, hue, k, max, maxHue, max_frac_widths_by_metric, max_val_by_metric, min, minHue, min_val_by_metric, rec, sampleReport, _i, _j, _k, _len, _len1, _len2, _ref, _ref1, _ref2, _results;
    max_frac_widths_by_metric = {};
    min_val_by_metric = {};
    max_val_by_metric = {};
    _ref = report.sample_reports;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      sampleReport = _ref[_i];
      calc_records_cell_contents(sampleReport.records, font);
      _ref1 = sampleReport.records;
      for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
        rec = _ref1[_j];
        if (!(rec.metric.name in max_frac_widths_by_metric)) {
          max_frac_widths_by_metric[rec.metric.name] = rec.frac_width;
        } else if (rec.frac_width > max_frac_widths_by_metric[rec.metric.name]) {
          max_frac_widths_by_metric[rec.metric.name] = rec.frac_width;
        }
        if (rec.num != null) {
          if (!(rec.metric.name in min_val_by_metric)) {
            min_val_by_metric[rec.metric.name] = rec.num;
          } else if (min_val_by_metric[rec.metric.name] > rec.num) {
            min_val_by_metric[rec.metric.name] = rec.num;
          }
          if (!(rec.metric.name in max_val_by_metric)) {
            max_val_by_metric[rec.metric.name] = rec.num;
          } else if (max_val_by_metric[rec.metric.name] < rec.num) {
            max_val_by_metric[rec.metric.name] = rec.num;
          }
        }
      }
    }
    _ref2 = report.sample_reports;
    _results = [];
    for (_k = 0, _len2 = _ref2.length; _k < _len2; _k++) {
      sampleReport = _ref2[_k];
      _results.push((function() {
        var _l, _len3, _ref3, _results1;
        _ref3 = sampleReport.records;
        _results1 = [];
        for (_l = 0, _len3 = _ref3.length; _l < _len3; _l++) {
          rec = _ref3[_l];
          if (rec.frac_width != null) {
            rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width;
            if (rec.right_shift !== 0) {
              a = 0;
            }
          }
          if (rec.num != null) {
            max = max_val_by_metric[rec.metric.name];
            min = min_val_by_metric[rec.metric.name];
            maxHue = GREEN_HUE;
            minHue = RED_HUE;
            if (rec.metric.quality === 'Less is better') {
              maxHue = RED_HUE;
              minHue = GREEN_HUE;
            }
            if (max === min) {

            } else {
              k = (maxHue - minHue) / (max - min);
              hue = Math.round(minHue + (rec.num - min) * k);
              _results1.push(rec.color = get_color(hue));
            }
          } else {
            _results1.push(void 0);
          }
        }
        return _results1;
      })());
    }
    return _results;
  };

  reporting.buildCommonRecords = function(common_records) {
    var rec, table, use_full_name, _i, _len;
    if (common_records) {
      calc_records_cell_contents(common_records, $('#report').css('font'));
      table = "<table cellspacing=\"0\" class=\"common_table\" id=\"common_table\">";
      for (_i = 0, _len = common_records.length; _i < _len; _i++) {
        rec = common_records[_i];
        table += "\n<tr><td> <span class='metric_name'>" + (get_metric_name_html(rec, use_full_name = true)) + ":</span> " + rec.cell_contents + " </td></tr>";
      }
      table += "\n</table>\n";
      return $('#report').append(table);
    }
  };

  reporting.buildTotalReport = function(report, columnOrder) {
    var padding, pos, rec, recNum, sampleLink, sampleName, sampleReport, table, _i, _j, _k, _len, _ref, _ref1, _ref2;
    if (report.name != null) {
      $('#report').append("<h3 class='table_name' style='margin: 0px 0 5px 0'>" + report.name + "</h3>");
    }
    calc_cell_contents(report, $('#report').css('font'));
    table = "<table cellspacing=\"0\" class=\"report_table " + (DRAGGABLE_COLUMNS ? 'draggable' : '') + " fix-align-char\" id=\"report_table_" + report.name + "\">";
    table += "\n<tr class=\"top_row_tr\">";
    table += "<td class=\"top_left_td left_column_td\"> <span>Sample</span> </td>";
    for (recNum = _i = 0, _ref = report.sample_reports[0].records.length; 0 <= _ref ? _i < _ref : _i > _ref; recNum = 0 <= _ref ? ++_i : --_i) {
      pos = columnOrder[recNum];
      rec = report.sample_reports[0].records[pos];
      table += "<td class='second_through_last_col_headers_td' position='" + pos + "'> {if DRAGGABLE_COLUMNS then '<span class=\'drag_handle\'><span class=\'drag_image\'></span></span>' else ''} <span class=\"metricName " + (DRAGGABLE_COLUMNS ? 'drag_handle' : '') + "\">" + (get_metric_name_html(rec)) + "</span> </td>";
    }
    _ref1 = report.sample_reports;
    for (_j = 0, _len = _ref1.length; _j < _len; _j++) {
      sampleReport = _ref1[_j];
      sampleName = sampleReport.sample.name;
      sampleLink = sampleReport.link;
      if (sampleName.length > 30) {
        sampleName = "<span title=\"" + sampleName + "\">" + (sampleName.trunc(80)) + "</span>";
      }
      table += "\n<tr> <td class=\"left_column_td\"> <a class=\"sample_name\" href=\"" + sampleLink + "\">" + sampleName + "</a> </td>";
      for (recNum = _k = 0, _ref2 = sampleReport.records.length; 0 <= _ref2 ? _k < _ref2 : _k > _ref2; recNum = 0 <= _ref2 ? ++_k : --_k) {
        pos = columnOrder[recNum];
        rec = sampleReport.records[pos];
        table += "<td metric=\"" + rec.metric.name + "\" style=\"" + CSS_PROP_TO_COLOR + ": " + rec.color + "\" class='number' quality=\"" + rec.metric.quality + "\"";
        if (rec.num != null) {
          table += ' number="' + rec.value + '">';
        }
        if (rec.right_shift != null) {
          padding = "margin-left: " + rec.right_shift + "px; margin-right: -" + rec.right_shift + "px;";
        } else {
          padding = "";
        }
        table += "<a style=\"" + padding + "\" " + (get_meta_tag_contents(rec)) + ">" + rec.cell_contents + " </a> </td>";
      }
      table += "</tr>";
    }
    table += "\n</table>\n";
    return $('#report').append(table);
  };

  set_legend = function() {
    var hue, legend, step, _i;
    legend = '<span>';
    step = 6;
    for (hue = _i = RED_HUE; step > 0 ? _i <= GREEN_HUE : _i >= GREEN_HUE; hue = _i += step) {
      legend += "<span style=\"" + CSS_PROP_TO_COLOR + ": " + (get_color(hue)) + "\">";
      switch (hue) {
        case RED_HUE:
          legend += 'w';
          break;
        case RED_HUE + step:
          legend += 'o';
          break;
        case RED_HUE + 2 * step:
          legend += 'r';
          break;
        case RED_HUE + 3 * step:
          legend += 's';
          break;
        case RED_HUE + 4 * step:
          legend += 't';
          break;
        case GREEN_HUE - 3 * step:
          legend += 'b';
          break;
        case GREEN_HUE - 2 * step:
          legend += 'e';
          break;
        case GREEN_HUE - step:
          legend += 's';
          break;
        case GREEN_HUE:
          legend += 't';
          break;
        default:
          legend += '.';
      }
      legend += "</span>";
    }
    legend += "</span>";
    return $('#report_legend').append(legend);
  };

  $.fn._splitDot_partTextWidth = function(html, font, part_type) {
    var frac_part, parts;
    parts = html.split('.');
    if (part_type === 'frac') {
      if (parts.length < 2) {
        return 0;
      } else {
        frac_part = '.' + parts[1];
      }
    } else if (part_type === 'int') {
      frac_part = parts[0];
    }
    if (!$.fn.fracPartTextWidth.fakeEl) {
      $.fn.fracPartTextWidth.fakeEl = $('<span>').hide().appendTo(document.body);
    }
    $.fn.fracPartTextWidth.fakeEl.html(frac_part);
    $.fn.fracPartTextWidth.fakeEl.css('font', font);
    return $.fn.fracPartTextWidth.fakeEl.width();
  };

  $.fn.fracPartTextWidth = function(html, font) {
    return $.fn._splitDot_partTextWidth(html, font, 'frac');
  };

  $.fn.intPartTextWidth = function(html, font) {
    return $.fn._splitDot_partTextWidth(html, font, 'int');
  };

  $.fn.textWidth = function(text, font) {
    if (!$.fn.textWidth.fakeEl) {
      $.fn.textWidth.fakeEl = $('<span>').hide().appendTo(document.body);
    }
    $.fn.textWidth.fakeEl.html(text);
    $.fn.textWidth.fakeEl.css('font', font);
    return $.fn.textWidth.fakeEl.width();
  };

  String.prototype.trunc = function(n) {
    var _ref;
    return this.substr(0, n - 1) + ((_ref = this.length > n) != null ? _ref : {
      '&hellip;': ''
    });
  };

}).call(this);

//# sourceMappingURL=build_total_report.map
