# help stuff

# in control panel director callback function
# to use, readd help_clicks as a parameter of the control panel director function
elif triggered_id == "help":
				return open_help(help_clicks)


# callbacks for after the panel director function

def open_help(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'help.n_clicks':
				return html.Div([
					dbc.Row([
						html.H4("Welcome to the help page", className="text-center", style={"margin-top": "20px", "color": "#102A5F"}),
						]),
					dbc.Row([
						html.Div([
							dbc.Button("pathways & reactions", id='pr_help',color="secondary", style={"margin-bottom": "10px", "margin-left": "70px",
								"margin-top": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							dbc.Button("element transfer network", id='etn_help',color="secondary", outline=True, style={
								"margin-left": "15px", "margin-bottom": "10px", "margin-top": "10px", "border-radius": "8px", "color": "#102A5F"}),
							]) 
						]),
					dbc.Row([
						html.Div([
							dbc.Button("gene delete", id='gd_help',color="secondary", outline=True, style={"margin-bottom": "10px",
								"margin-left": "32px", "border-radius": "8px", "color": "#102A5F"},),
							dbc.Button("path & compound search", id='pcs_help',color="secondary", style={"margin-left": "15px", 
								"margin-bottom": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							dbc.Button("flux analysis", id='fa_help',color="secondary", style={"margin-left": "15px", 
								"margin-bottom": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							])
						]),
					dbc.Row([html.Div(id="help_area", className="help", children=[])]),
					])

		@_app.callback(
			Output("help_area", "children"),
			Input("pr_help", "n_clicks"),
			Input("etn_help", "n_clicks"),
			Input("gd_help", "n_clicks"),
			Input("pcs_help", "n_clicks"),
			Input("fa_help", "n_clicks"),
		)
		def det_help(pr_clicks, etn_clicks, gd_clicks, pcs_clicks, fa_clicks):
			triggered_id = ctx.triggered_id
			if triggered_id == "pr_help":
				return open_pr_help(pr_clicks)
			elif triggered_id == "etn_help":
				return open_etn_help(etn_clicks)
			elif triggered_id == "gd_help":
				return open_gd_help(gd_clicks)
			elif triggered_id == "pcs_help":
				return open_pcs_help(pcs_clicks)
			elif triggered_id == "fa_help":
				return open_fa_help(fa_clicks)
		def open_pr_help(n_clicks):
			return html.Div([
				dbc.Row([
					html.Div([
						html.H4("Pathways", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the dropdown to display pathways in the metabolic model. Multiple can be selected.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				dbc.Row([
					html.Div([
						html.H4("Reactions", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the dropdown to display reactions in the metabolic model. Multiple can be selected.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				])
		def open_etn_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H4("Element Transfer Networks", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the dropdown to choose a chemical element for the network (default Carbon).",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				])
		def open_gd_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H4("Gene Delete", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("""Use the dropdown to delete genes from the model. If all genes associated with a reaction are
							deleted, the reaction will be removed from analysis.""",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				])
		def open_pcs_help(nclicks):
			return html.Div([
				dbc.Row([
					html.Div([
						html.H4("Path Search", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the path search to conduct a bidirectional breadth first search to show the shortest path between two compounds.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				dbc.Row([
					html.Div([
						html.H4("Compound Search", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the compound search to show all reactions containing this compound. Select ALL for pathways to see all reactions.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				])
		def open_fa_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H4("Flux Analysis", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H5("Use the dropdown to choose a reaction to optimize for flux balance analysis. Flux is visualized in blue.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F",
							"text-transform": "none"}),
						]),
					]),
				])


"""dbc.Button("help", id="help", color="primary", style={"margin-left": "4px", "margin-bottom": "10px",
											"background-color": "#102A5F", "border-radius": "8px", "margin-top": "10px"},),"""