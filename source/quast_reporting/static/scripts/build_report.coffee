showPlotWithInfo = (info) ->
  newSeries = []
  newColors = []

  $('#legend-placeholder').find 'input:checked' .each ->
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


recoverOrderFromCookies = ->
  return null unless navigator.cookieEnabled

  orderString = readCookie "order"
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
    result = null;

  return result;


totalReportData =
  report: null
  date: null
  order: null


reporting.buildReport = ->
  unless (totalReportData = readJson('total-report'))
    console.log "Error: cannot read #total-report-json"
    return 1

  report = totalReportData.report
  return 1 if report is []

  date = totalReportData.date
  order = totalReportData.order

  columnNames = (metric.name for metric in report[0].metrics)
  columnOrder = recoverOrderFromCookies() or order or [0...columnNames.length]

  reporting.buildTotalReport(report, columnOrder, date)

  return 0

