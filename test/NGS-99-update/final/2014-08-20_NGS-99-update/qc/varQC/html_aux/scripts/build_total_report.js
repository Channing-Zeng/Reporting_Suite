//@ sourceMappingURL=build_total_report.map
// Generated by CoffeeScript 1.6.1
(function() {
  var BLUE_HUE, CSS_PROP_TO_COLOR, DRAGGABLE_COLUMNS, GREEN_HSL, GREEN_HUE, RED_HUE, calc_cell_contents, calc_records_cell_contents, check_all_values_equal, get_color, get_meta_tag_contents, get_metric_name_html, mean, median, metric, record, sampleReport, set_legend,
    __hasProp = {}.hasOwnProperty;

  sampleReport = {
    sample: {
      name: '',
      display_name: '',
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
    html_fpath: '',
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
    type: null,
    all_values_equal: false
  };

  DRAGGABLE_COLUMNS = false;

  RED_HUE = 0;

  GREEN_HUE = 120;

  BLUE_HUE = 240;

  GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)';

  CSS_PROP_TO_COLOR = 'background-color';

  get_color = function(hue, lightness) {
    lightness = lightness != null ? lightness : 92;
    return 'hsl(' + hue + ', 80%, ' + lightness + '%)';
  };

  check_all_values_equal = function(vals) {
    var first_val, val, _i, _len;
    first_val = null;
    for (_i = 0, _len = vals.length; _i < _len; _i++) {
      val = vals[_i];
      if (first_val != null) {
        if (val !== first_val) {
          return false;
        }
      } else {
        first_val = val;
      }
    }
    return true;
  };

  get_meta_tag_contents = function(rec) {
    var a, db, dbs, k, meta, meta_table, novelty, short_table, val, val_by_db, values, _i, _len;
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
          val_by_db = meta[novelty];
          dbs = (function() {
            var _results;
            _results = [];
            for (db in val_by_db) {
              val = val_by_db[db];
              if (db !== 'average') {
                _results.push(db);
              }
            }
            return _results;
          })();
          dbs.push('average');
          break;
        }
        short_table = true;
        for (novelty in meta) {
          val_by_db = meta[novelty];
          if (!check_all_values_equal((function() {
            var _results;
            _results = [];
            for (db in val_by_db) {
              val = val_by_db[db];
              if (db !== 'average') {
                _results.push(val);
              }
            }
            return _results;
          })())) {
            short_table = false;
          }
        }
        if (short_table) {
          meta_table += '<tr><td></td>';
          for (novelty in meta) {
            val_by_db = meta[novelty];
            if (novelty !== 'all') {
              meta_table += "<td>" + (toPrettyString(val_by_db[dbs[0]], rec.metric.unit)) + "</td>";
            }
          }
          meta_table += '</tr>\n';
        } else {
          for (_i = 0, _len = dbs.length; _i < _len; _i++) {
            db = dbs[_i];
            meta_table += "<tr><td>" + db + "</td>";
            for (novelty in meta) {
              val_by_db = meta[novelty];
              if (novelty !== 'all') {
                meta_table += "<td>" + (toPrettyString(val_by_db[db], rec.metric.unit)) + "</td>";
              }
            }
            meta_table += '</tr>\n';
          }
        }
        meta_table += '</table>\n';
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"" + meta_table + "\"";
      }
    } else {
      return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\"";
    }
  };

  get_metric_name_html = function(metric, use_full_name) {
    var description, metricName;
    if (use_full_name == null) {
      use_full_name = false;
    }
    if (metric.short_name && !use_full_name) {
      metricName = metric.short_name;
      description = metric.description || metric.name;
      return "<a class=\"tooltip-link\" rel=\"tooltip\" title=\"" + description + "\">" + metricName + "</a>";
    } else {
      return metric.name;
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

  mean = function(a, b) {
    return (a + b) / 2;
  };

  calc_cell_contents = function(report, section, font) {
    var a, l, max_bri, max_frac_widths_by_metric, min_bri, rec, vals, _i, _j, _k, _l, _len, _len1, _len2, _len3, _len4, _m, _ref, _ref1, _ref2, _ref3, _ref4;
    max_frac_widths_by_metric = {};
    _ref = report.sample_reports;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      sampleReport = _ref[_i];
      calc_records_cell_contents(sampleReport.records, font);
      _ref1 = sampleReport.records;
      for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
        rec = _ref1[_j];
        if (!(rec.metric.name in section.metrics_by_name)) {
          continue;
        }
        if (!(rec.metric.name in max_frac_widths_by_metric)) {
          max_frac_widths_by_metric[rec.metric.name] = rec.frac_width;
        } else if (rec.frac_width > max_frac_widths_by_metric[rec.metric.name]) {
          max_frac_widths_by_metric[rec.metric.name] = rec.frac_width;
        }
        if (rec.num != null) {
          if (rec.metric.values == null) {
            rec.metric.values = [];
          }
          rec.metric.values.push(rec.num);
        }
      }
    }
    _ref2 = section.metrics;
    for (_k = 0, _len2 = _ref2.length; _k < _len2; _k++) {
      metric = _ref2[_k];
      if (!(metric.values != null)) {
        continue;
      }
      vals = metric.values.slice().sort(function(a, b) {
        return a - b;
      });
      l = vals.length;
      metric.min = vals[0];
      metric.max = vals[vals.length - 1];
      metric.med = l % 2 !== 0 ? vals[(l - 1) / 2] : mean(vals[l / 2], vals[(l / 2) - 1]);
      metric.q1 = vals[Math.floor((l - 1) / 4)];
      metric.q3 = vals[Math.floor((l - 1) * 3 / 4)];
      metric.d = metric.q3 - metric.q1;
    }
    _ref3 = report.sample_reports;
    for (_l = 0, _len3 = _ref3.length; _l < _len3; _l++) {
      sampleReport = _ref3[_l];
      _ref4 = sampleReport.records;
      for (_m = 0, _len4 = _ref4.length; _m < _len4; _m++) {
        rec = _ref4[_m];
        if (!(rec.metric.name in section.metrics_by_name)) {
          continue;
        }
        if (rec.frac_width != null) {
          rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width;
          if (rec.right_shift !== 0) {
            a = 0;
          }
        }
        metric = rec.metric;
        if (rec.num != null) {
          max_bri = 0;
          min_bri = 100;
          if (metric.min === metric.max) {
            metric.all_values_equal = true;
          } else {
            metric.all_values_equal = false;
            if (rec.num < metric.q1 - 3 * metric.d) {
              rec.color = get_color(BLUE_HUE, 30);
            } else if (rec.num < metric.q1 - 1.5 * metric.d) {
              rec.color = get_color(BLUE_HUE, 60);
            } else if (rec.num > metric.q3 + 1.5 * metric.d) {
              rec.color = '#88FFFF';
            } else if (rec.num > metric.q3 + 3 * metric.d) {
              rec.color = '#55FFFF';
            }
          }
        }
      }
    }
    return report;
  };

  median = function(x) {
    var sorted;
    if (x.length === 0) {
      return null;
    }
    sorted = x.slice().sort(function(a, b) {
      return a - b;
    });
    if (sorted.length % 2 === 1) {
      return sorted[(sorted.length - 1) / 2];
    } else {
      return (sorted[(sorted.length / 2) - 1] + sorted[sorted.length / 2]) / 2;
    }
  };

  reporting.buildTotalReport = function(report, section, columnOrder) {
    var colNum, i, line_caption, padding, pos, r, rec, sort_by, table, _i, _j, _k, _l, _len, _len1, _ref, _ref1, _ref2, _ref3;
    if (section.title != null) {
      $('#report').append("<h3 class='table_name' style='margin: 0px 0 5px 0'>" + section.title + "</h3>");
    }
    calc_cell_contents(report, section, $('#report').css('font'));
    table = "<table cellspacing=\"0\"                    class=\"report_table tableSorter " + (DRAGGABLE_COLUMNS ? 'draggable' : '') + " fix-align-char\"                    id=\"report_table_" + section.name + "\">";
    table += "\n<tr class=\"top_row_tr\">";
    table += "<th class=\"top_left_td left_column_td\" data-sortBy='numeric'>                    <span>Sample</span>              </th>";
    for (colNum = _i = 0, _ref = section.metrics.length; 0 <= _ref ? _i < _ref : _i > _ref; colNum = 0 <= _ref ? ++_i : --_i) {
      pos = columnOrder[colNum];
      metric = section.metrics[pos];
      sort_by = metric.all_values_equal ? 'nosort' : 'numeric';
      table += "<th class='second_through_last_col_headers_td' data-sortBy=" + sort_by + " position='" + pos + "'>             <span class=\'metricName " + (DRAGGABLE_COLUMNS ? 'drag_handle' : '') + "\'>" + (get_metric_name_html(metric)) + "</span>        </th>";
    }
    i = 0;
    _ref1 = report.sample_reports;
    for (_j = 0, _len = _ref1.length; _j < _len; _j++) {
      sampleReport = _ref1[_j];
      line_caption = sampleReport.display_name;
      if (line_caption.length > 30) {
        line_caption = "<span title=\"" + line_caption + "\">" + (line_caption.trunc(80)) + "</span>";
      }
      table += "\n<tr>            <td class=\"left_column_td\" data-sortAs=" + (report.sample_reports.length - i) + ">";
      if (report.sample_reports.length === 1) {
        table += "<span class=\"sample_name\">" + line_caption + "</span>";
      } else {
        if (sampleReport.html_fpath != null) {
          table += "<a class=\"sample_name\" href=\"" + sampleReport.html_fpath + "\">" + line_caption + "</a>";
        } else {
          table += "<span class=\"sample_name\"\">" + line_caption + "</span>";
        }
      }
      table += "</td>";
      for (colNum = _k = 0, _ref2 = section.metrics.length; 0 <= _ref2 ? _k < _ref2 : _k > _ref2; colNum = 0 <= _ref2 ? ++_k : --_k) {
        pos = columnOrder[colNum];
        metric = section.metrics[pos];
        rec = null;
        _ref3 = sampleReport.records;
        for (_l = 0, _len1 = _ref3.length; _l < _len1; _l++) {
          r = _ref3[_l];
          if (r.metric.name === metric.name) {
            rec = r;
            break;
          }
        }
        if (rec == null) {
          table += "<td></td>";
          continue;
        }
        table += "<td metric=\"" + metric.name + "\"                          style=\"" + CSS_PROP_TO_COLOR + ": " + rec.color + "\"                          class='number'                          quality=\"" + metric.quality + "\"";
        if (rec.num != null) {
          table += " number=\"" + rec.value + "\" data-sortAs=" + rec.value + ">";
        } else {
          table += ">";
        }
        if (rec.right_shift != null) {
          padding = "margin-left: " + rec.right_shift + "px; margin-right: -" + rec.right_shift + "px;";
        } else {
          padding = "";
        }
        table += "<a style=\"" + padding + "\"                          " + (get_meta_tag_contents(rec)) + ">" + rec.cell_contents + "                      </a>                    </td>";
      }
      table += "</tr>";
      i += 1;
    }
    table += "\n</table>\n";
    return $('#report').append(table);
  };

  reporting.buildCommonRecords = function(common_records) {
    var rec, table, use_full_name, _i, _len;
    if (common_records.length === 0) {
      return;
    }
    calc_records_cell_contents(common_records, $('#report').css('font'));
    table = "<table cellspacing=\"0\" class=\"common_table\" id=\"common_table\">";
    for (_i = 0, _len = common_records.length; _i < _len; _i++) {
      rec = common_records[_i];
      table += "\n<tr><td>                <span class='metric_name'>" + (get_metric_name_html(rec.metric, use_full_name = true)) + ":</span>                " + rec.cell_contents + "              </td></tr>";
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
