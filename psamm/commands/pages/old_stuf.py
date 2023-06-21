# Define the main body of the app
    body_layout = dbc.Container(
        [dbc.Row([
            dbc.Row([
                html.Div(id='modelling', className='modelling', children=[
                    dcc.Tabs(id="parent-tabs", children=[
                        dcc.Tab(id="sim", label='Simulate', style={
                            'width': '30%', 'font-size': '150%', 'height': '50%'}, children=[
                            dbc.Row([
                                dbc.Col([
                                    html.Div(id='model-control-tabs',
                                             className='control-tabs',
                                             style={'width': '69%', 'font-size': '50%', 
                                             'height': '50%'}, 
                                              children=[
                                                 simulate_tabs 
                                             ]), ]),
                                dbc.Col([
                                    s_n,
                                    dbc.Row([
                                        dbc.Alert(
                                            id="node-data",
                                            children="Click on a node to "
                                                     "see its details",
                                            color="secondary",
                                        ),
                                    ]), ]),
                                dbc.Row([dcc.Markdown(
                                    '''
                                    \\* Data analysis carried out for
                                    demonstration of data visualisation
                                    purposes only.
                                    ''')], style={"fontSize": 11,
                                                  "color": "gray"},), ],
                                    style={"marginTop": 20},)]),
                        dcc.Tab(label='Curate', id='cur', style={
                            'width': '100%', 'font-size': '150%', 'height': '50%'},
                                children=[
                                    html.Div(id='curate-control-tabs',
                                             className='curate-tabs',
                                             children=[dcc.Tabs(
                                                id="curate-subtabs",
                                                value="what is",
                                                children=[add_tab,
                                                          e_t,
                                                          chargecheck,
                                                          fc, ]),

                                                       ]),
                                    dbc.Row([
                                        dbc.Alert(
                                            id="save_confirmation",
                                            children="Press Save Model "
                                            "to save the model",
                                            color="secondary",
                                        ),
                                        dbc.Row([
                                            dbc.Col([
                                                html.Button("Save Model",
                                                    id="btn_save_model"),
                                            ]), ]), ],
                                        style={"fontSize": 16,
                                               "color": "gray"},),
                                ]),
                    ]),
                ]), ]),
        ]),
        ], 
        # take up the whole screen
        fluid = True, style={"height": "100vh"} 
        )