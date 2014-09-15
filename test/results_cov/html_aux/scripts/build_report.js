// Generated by CoffeeScript 1.7.1
(function() {
  var extend, merge, metric, preprocessReport, readJson, recoverOrderFromCookies, report, section, showPlotWithInfo, totalReportData;

  showPlotWithInfo = function(info) {
    var newColors, newSeries;
    newSeries = [];
    newColors = [];
    return $('#legend-placeholder').find('input:checked'.each(function() {
      var i, number, series, _i, _ref;
      number = $(this).attr('name');
      if (number && info.series && info.series.length > 0) {
        for (i = _i = i, _ref = info.series.length; i <= _ref ? _i < _ref : _i > _ref; i = i <= _ref ? ++_i : --_i) {
          series = info.series[i];
          if (series.number !== number) {
            break;
          }
        }
        if (i <= info.series.length) {
          newSeries.push(series);
          newColors.push(series.color);
        } else {
          console.log('no series with number ' + number);
        }
      }
      if (newSeries.length === 0) {
        newSeries.push({
          data: []
        });
        newColors.push('#FFF');
      }
      return info.showWithData(newSeries, newColors);
    }));
  };

  recoverOrderFromCookies = function(report_name) {
    var columnOrder, fail, orderString, val, _i, _len, _ref;
    if (!navigator.cookieEnabled) {
      return null;
    }
    orderString = readCookie(report_name + '_order');
    if (!orderString) {
      return null;
    }
    columnOrder = [];
    fail = false;
    _ref = orderString.split(' ');
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      val = _ref[_i];
      val = parseInt(val);
      if (isNaN(val)) {
        fail = true;
      } else {
        columnOrder.push(val);
      }
    }
    if (fail) {
      return null;
    }
    return columnOrder;
  };

  readJson = function(what) {
    result;
    var e, result;
    try {
      result = JSON.parse($('#' + what + '-json').html());
    } catch (_error) {
      e = _error;
      result = null;
    }
    return result;
  };

  totalReportData = {
    date: null,
    report: null
  };

  report = {
    name: '',
    order: null,
    sample_reports: [],
    metric_storage: {
      common_for_all_samples_section: {
        name: '',
        metrics: []
      },
      sections: []
    }
  };

  section = {
    name: '',
    title: '',
    metrics: [],
    metrics_by_name: {}
  };

  metric = {
    name: '',
    short_name: '',
    description: '',
    quality: '',
    common: true,
    unit: ''
  };

  extend = function(object, properties) {
    var key, val;
    for (key in properties) {
      val = properties[key];
      object[key] = val;
    }
    return object;
  };

  merge = function(options, overrides) {
    return extend(extend({}, options), overrides);
  };

  preprocessReport = function(report) {
    var all_metrics_by_name, rec, s, sample_report, _i, _j, _k, _len, _len1, _len2, _ref, _ref1, _ref2;
    all_metrics_by_name = {};
    extend(all_metrics_by_name, report.metric_storage.common_for_all_samples_section.metrics_by_name);
    _ref = report.metric_storage.sections;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      s = _ref[_i];
      extend(all_metrics_by_name, s.metrics_by_name);
    }
    _ref1 = report.sample_reports;
    for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
      sample_report = _ref1[_j];
      sample_report.metric_storage = report.metric_storage;
      _ref2 = sample_report.records;
      for (_k = 0, _len2 = _ref2.length; _k < _len2; _k++) {
        rec = _ref2[_k];
        rec.metric = all_metrics_by_name[rec.metric.name];
      }
    }
    return report;
  };

  reporting.buildReport = function() {
    var columnNames, columnOrder, common_metrics_by_name, general_records, m, plot, plots_html, rec, sample_report, sample_reports, _i, _j, _k, _l, _len, _len1, _len2, _ref, _ref1, _ref2, _results;
    if (!(totalReportData = readJson('total-report'))) {
      console.log("Error: cannot read #total-report-json");
      return 1;
    }
    report = preprocessReport(totalReportData.report);
    $('#report_date').html('<p>' + totalReportData.date + '</p>');
    common_metrics_by_name = report.metric_storage.common_for_all_samples_section.metrics_by_name;
    general_records = (function() {
      var _i, _len, _ref, _results;
      _ref = report.sample_reports[0].records;
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        rec = _ref[_i];
        if (rec.metric.name in common_metrics_by_name) {
          _results.push(rec);
        }
      }
      return _results;
    })();
    reporting.buildCommonRecords(general_records);
    sample_reports = report.sample_reports;
    _ref = report.metric_storage.sections;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      section = _ref[_i];
      columnNames = (function() {
        var _j, _len1, _ref1, _results;
        _ref1 = section.metrics;
        _results = [];
        for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
          m = _ref1[_j];
          _results.push(m.name);
        }
        return _results;
      })();
      columnOrder = (recoverOrderFromCookies(section.name)) || report.order || (function() {
        _results = [];
        for (var _j = 0, _ref1 = columnNames.length; 0 <= _ref1 ? _j < _ref1 : _j > _ref1; 0 <= _ref1 ? _j++ : _j--){ _results.push(_j); }
        return _results;
      }).apply(this);
      reporting.buildTotalReport(report, section, columnOrder);
      plots_html = "";
      for (_k = 0, _len1 = sample_reports.length; _k < _len1; _k++) {
        sample_report = sample_reports[_k];
        _ref2 = sample_report.plots;
        for (_l = 0, _len2 = _ref2.length; _l < _len2; _l++) {
          plot = _ref2[_l];
          plots_html += "<img src=\"" + plot + "\"/>";
        }
      }
      $('#plot').html(plots_html);
    }
    return 0;
  };

}).call(this);

//# sourceMappingURL=build_report.map
