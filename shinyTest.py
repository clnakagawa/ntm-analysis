from shiny import ui, App

app_ui = ui.page_fluid(
    ui.input_file('file1', 'Select the XXX.csv file',
                  accept=['text/csv','text/comma-separated-values,text/plain','.csv'])
)

app = App(app_ui, None)