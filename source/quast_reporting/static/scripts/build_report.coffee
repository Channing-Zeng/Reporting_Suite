showPlotWithInfo = (info) ->
    newSeries = []
    newColors = []

    $('#legend-placeholder').find 'input:checked'.each ->
        number = $(this).attr 'name'
        if number and info.series && info.series.length > 0
            for i in [i...info.series.length]
                series = info.series[i]
                break if series.number != number

            if i <= info.series.length
                newSeries.push(series)
                newColors.push(series.color)
            else
                console.log('no series with number ' + number)

        if newSeries.length == 0
            newSeries.push
                data: []

            newColors.push '#FFF'

        info.showWithData(newSeries, newColors)


recoverOrderFromCookies = (report_name) ->
    return null unless navigator.cookieEnabled

    orderString = readCookie report_name + '_order'
    return null unless orderString

    columnOrder = []
    fail = false

    for val in orderString.split(' ')
        val = parseInt val
        if isNaN val
            fail = true
        else
            columnOrder.push val

    return null if fail

    return columnOrder


readJson = (what) ->
    result
    try
        result = JSON.parse $('#' + what + '-json').html()
    catch e
        result = null

    return result


totalReportData =
    date: null
    report: null

report =
    name: ''
    order: null
    sample_reports: []
    metric_storage:
        common_for_all_samples_section:
            name: ''
            metrics: []
        sections: []

section =
    name: ''
    title: ''
    metrics: []
    metrics_by_name: {}

metric =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    common: true
    unit: ''


extend = (object, properties) ->
    for key, val of properties
        object[key] = val
    object


merge = (options, overrides) ->
    extend (extend {}, options), overrides


preprocessReport = (report) ->
    all_metrics_by_name = {}
    extend all_metrics_by_name, report.metric_storage.common_for_all_samples_section.metrics_by_name
    for s in report.metric_storage.sections
        extend all_metrics_by_name, s.metrics_by_name

    for sample_report in report.sample_reports
        sample_report.metric_storage = report.metric_storage
        for rec in sample_report.records
            rec.metric = all_metrics_by_name[rec.metric.name]

    return report


reporting.buildReport = ->
    unless (totalReportData = readJson 'total-report')
        console.log "Error: cannot read #total-report-json"
        return 1

    report = preprocessReport totalReportData.report

    $('#report_date').html '<p>' + totalReportData.date + '</p>'

    common_metrics_by_name = report.metric_storage.common_for_all_samples_section.metrics_by_name
    general_records = (rec for rec in report.sample_reports[0].records when rec.metric.name of common_metrics_by_name)
    reporting.buildCommonRecords general_records

    sample_reports = report.sample_reports
    for section in report.metric_storage.sections
        columnNames = (r.metric.name for r in sample_reports[0].records when r.metric.name of section.metrics_by_name)
        columnOrder = (recoverOrderFromCookies section.name) or report.order or [0...columnNames.length]

        reporting.buildTotalReport report, section, columnOrder
        plots_html = ""
        for sample_report in sample_reports
            for plot in sample_report.plots
                plots_html += "<img src=\"#{plot}\"/>"
        $('#plot').html plots_html

    return 0