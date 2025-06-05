for gene in genes_detected:
    gene_name_key = gene["name"].strip('"')
    # Usar una lista ordenada para mantener el orden en la UI y al guardar
    proteins_for_gene = sorted(list({p.strip('"') for p in gene["proteins"]}))

    # VBox principal para la UI de este gen
    gene_box_content = []
    
    # 1. Encabezado del gen (nombre y proteínas que produce)
    gene_info_html = f"<b>Gene: {gene_name_key}</b> (Produces: {', '.join(proteins_for_gene)})"
    gene_box_content.append(widgets.HTML(value=gene_info_html))

    # Diccionario para almacenar los widgets de parámetros de cada proteína
    protein_params_widgets = {}
    
    # 2. Lógica para los parámetros de tiempo (depende del número de proteínas)
    if len(proteins_for_gene) == 1:
        # --- CASO 1: UNA SOLA PROTEÍNA (Diseño original) ---
        protein_id = proteins_for_gene[0]
        gene_box_content.append(widgets.Label("Timing (Activation/Degradation, minutes):"))
        
        act_time_w = widgets.FloatText(description="Act. Time:", value=15.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
        act_var_w = widgets.FloatText(description="Act. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})
        deg_time_w = widgets.FloatText(description="Deg. Time:", value=20.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
        deg_var_w = widgets.FloatText(description="Deg. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})

        protein_params_widgets[protein_id] = {
            "act_time": act_time_w, "act_var": act_var_w,
            "deg_time": deg_time_w, "deg_var": deg_var_w
        }
        protein_hbox = widgets.HBox([act_time_w, create_spacer(), act_var_w, create_spacer(), deg_time_w, create_spacer(), deg_var_w], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
        gene_box_content.append(protein_hbox)

    else:
        # --- CASO 2: MÚLTIPLES PROTEÍNAS ---
        gene_box_content.append(widgets.Label("Timing (Activation/Degradation, minutes):"))
        for protein_id in proteins_for_gene:
            # Encabezado para cada proteína
            protein_header = widgets.HTML(value=f"<b>Protein: {protein_id}</b>", layout=widgets.Layout(margin='5px 0 2px 0'))
            gene_box_content.append(protein_header)

            # Widgets con descripciones genéricas, ya que el encabezado especifica la proteína
            act_time_w = widgets.FloatText(description="Act. Time:", value=15.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
            act_var_w = widgets.FloatText(description="Act. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})
            deg_time_w = widgets.FloatText(description="Deg. Time:", value=20.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
            deg_var_w = widgets.FloatText(description="Deg. Var.:", value=1.0, layout=content_widget_layout_small, style={'description_width': 'initial'})

            protein_params_widgets[protein_id] = {
                "act_time": act_time_w, "act_var": act_var_w,
                "deg_time": deg_time_w, "deg_var": deg_var_w
            }
            protein_hbox = widgets.HBox([act_time_w, create_spacer(), act_var_w, create_spacer(), deg_time_w, create_spacer(), deg_var_w], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
            gene_box_content.append(protein_hbox)

    # 3. Parámetros de ruido (comunes para todo el gen)
    gene_box_content.append(widgets.Label("Noise Parameters:"))
    noise_widgets = {
      "toOn": widgets.FloatSlider(description="P(noise ON):", min=0.0, max=1.0, step=0.01, value=0.0, layout=content_widget_layout_medium, style={'description_width': 'initial'}),
      "toOff": widgets.FloatSlider(description="P(noise Off):", min=0.0, max=1.0, step=0.01, value=0.0, layout=content_widget_layout_medium, style={'description_width': 'initial'}),
      "noise_time": widgets.FloatText(description="Noise Act. Time:", value=100.0, layout=content_widget_layout_medium, style={'description_width': 'initial'})
    }
    noise_params_hbox = widgets.HBox([
      noise_widgets["noise_time"], create_spacer(),
      noise_widgets["toOn"], create_spacer(),
      noise_widgets["toOff"]
    ], layout=widgets.Layout(width='100%', justify_content='flex-start', align_items='center'))
    gene_box_content.append(noise_params_hbox)

    # 4. Guardar los widgets y ensamblar la UI para este gen
    gene_widgets[gene_name_key] = {
        "protein_params": protein_params_widgets,
        "proteins_order": proteins_for_gene,
        "noise_params": noise_widgets
    }
    
    gene_box = widgets.VBox(gene_box_content, layout=widgets.Layout(border='1px dashed lightgray', margin='5px 0', padding='5px', width='100%'))
    all_genes_boxes.append(gene_box)

  genes_section_box = widgets.VBox([widgets.HTML(value="<h2>Gene Parameters</h2>")] + all_genes_boxes,
                                   layout=widgets.Layout(border='1px solid lightgray', margin='10px 0', padding='10px', width='100%'))
