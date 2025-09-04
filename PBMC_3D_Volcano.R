"R Script for the visualisation of regulated proteins compared to baseline separately in PRO and CTR in PBMCs.
Output: 3D Volcanos"


# === Load required packages ===
install.packages("plotly")

library(data.table)
library(reshape2)
library(plotly)
library(htmlwidgets)

# === Load dataset ===
file <- "/Users/majaewen/Documents/Uni/6_Semester/Daten/Volcanos/T_Test_Profis.txt" #Enter T_Test_Controls.txt to generate the Volcanos for CTR
prot_toplist <- fread(file)
setDT(prot_toplist)


cols_to_fix <- c(
  "Welch's T-test Difference T1_T0", 
  "Welch's T-test Difference T2_T0", 
  "Welch's T-test Difference T3_T0", 
  "Welch's T-test Difference T4_T0",
  "-Log Welch's T-test p-value T1_T0",
  "-Log Welch's T-test p-value T2_T0",
  "-Log Welch's T-test p-value T3_T0",
  "-Log Welch's T-test p-value T4_T0",
  "Welch's T-test q-value T1_T0",
  "Welch's T-test q-value T2_T0",
  "Welch's T-test q-value T3_T0",
  "Welch's T-test q-value T4_T0"
)


prot_toplist[, (cols_to_fix) := lapply(.SD, function(x) {
  as.numeric(gsub(",", ".", x))
}), .SDcols = cols_to_fix]



setDT(prot_toplist)

# === Reshape -log10(p-values) ===
p_value_long <- melt(prot_toplist, 
                     id.vars = c("PG.ProteinGroups", "PG.Genes"), 
                     measure.vars = c("-Log Welch's T-test p-value T1_T0", 
                                      "-Log Welch's T-test p-value T2_T0", 
                                      "-Log Welch's T-test p-value T3_T0", 
                                      "-Log Welch's T-test p-value T4_T0"),
                     variable.name = "Comparison", 
                     value.name = "logP")

# === Reshape adjusted p-values (q-values) ===
adj_p_value_long <- melt(prot_toplist, 
                         id.vars = c("PG.ProteinGroups", "PG.Genes"), 
                         measure.vars = c("Welch's T-test q-value T1_T0", 
                                          "Welch's T-test q-value T2_T0", 
                                          "Welch's T-test q-value T3_T0", 
                                          "Welch's T-test q-value T4_T0"),
                         variable.name = "Comparison", 
                         value.name = "adj.P.Val")

# === Reshape log fold changes ===
estimate_long <- melt(prot_toplist, 
                      id.vars = c("PG.ProteinGroups", "PG.Genes"), 
                      measure.vars = c("Welch's T-test Difference T1_T0", 
                                       "Welch's T-test Difference T2_T0", 
                                       "Welch's T-test Difference T3_T0", 
                                       "Welch's T-test Difference T4_T0"),
                      variable.name = "Comparison", 
                      value.name = "logFC")

# === Clean Comparison column ===
setDT(p_value_long)
setDT(adj_p_value_long)
setDT(estimate_long)

p_value_long[, Comparison := gsub("-Log Welch's T-test p-value ", "", Comparison)]
adj_p_value_long[, Comparison := gsub("Welch's T-test q-value ", "", Comparison)]
estimate_long[, Comparison := gsub("Welch's T-test Difference ", "", Comparison)]

# === Merge all reshaped data ===
prot_toplist_long <- merge(p_value_long, adj_p_value_long, 
                           by = c("PG.ProteinGroups", "PG.Genes", "Comparison"))
prot_toplist_long <- merge(prot_toplist_long, estimate_long, 
                           by = c("PG.ProteinGroups", "PG.Genes", "Comparison"))

# === Rename for consistency ===
setnames(prot_toplist_long, c("PG.ProteinGroups", "PG.Genes"), c("ID", "Gene"))

# === Ensure all numeric fields ===
prot_toplist_long[, logP := suppressWarnings(as.numeric(logP))]
prot_toplist_long[, logFC := suppressWarnings(as.numeric(logFC))]
prot_toplist_long[, adj.P.Val := suppressWarnings(as.numeric(adj.P.Val))]

# === Compute p-value and neg_log10(p-value) ===
prot_toplist_long[, P.Value := 10^(-logP)]
prot_toplist_long[, neg_log10_P.Value := logP]

# === Assign time point (exact match) ===
prot_toplist_long[, timepoint := fifelse(Comparison == "T1_T0", 1,
                                         fifelse(Comparison == "T2_T0", 2,
                                                 fifelse(Comparison == "T3_T0", 3,
                                                         fifelse(Comparison == "T4_T0", 4, NA_integer_))))]

# === Label for hover ===
prot_toplist_long[, lbl := paste0(
  "Gene: ", Gene,
  "<br>Protein: ", ID,
  "<br>logFC: ", sprintf("%.2e", logFC),
  "<br>p-value: ", sprintf("%.2e", P.Value),
  "<br>adj. p-value: ", sprintf("%.2e", adj.P.Val)
)]

# === Significance and color coding ===
prot_toplist_long[, sig := fcase(
  adj.P.Val < 0.05 & logFC > 0, "sig_up",
  adj.P.Val < 0.05 & logFC < 0, "sig_down",
  default = "not_sig"
)]

prot_toplist_long[, cols := fcase(
  sig == "sig_up", "#B22222",
  sig == "sig_down", "#00008B",
  default = "grey"
)]

# === Filter out rows with NA or invalid values ===
prot_toplist_filtered <- prot_toplist_long[
  !is.na(timepoint) & is.finite(logFC) & is.finite(neg_log10_P.Value)
]

# === Create 3D plot ===
phos3D_J <- plot_ly(data = prot_toplist_filtered) %>%
  add_trace(
    x = ~ timepoint, 
    y = ~ logFC,
    z = ~ neg_log10_P.Value,
    text = ~ lbl,
    hoverinfo = 'text',
    type = 'scatter3d', 
    mode = 'markers',
    marker = list(
      color = ~ cols,
      size = 3,
      line = list(color = 'whitesmoke', width = 0.1)
    )
  ) %>%
  layout(
    annotations = list(
      list(
        text = "<b>PRO - PBMCs</b>",
        x = 0.5,
        y = 0.85,  
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
       font = list(size = 30, color = "black")
      )
    ),
    margin = list(t = 30),
    scene = list(
      xaxis = list(
        tickmode = "array",
        title = list(text = "Timepoint", font = list(size = 18, color = "black")),
        tickvals = c(1, 2, 3, 4),
        ticktext = c("1", "2", "3", "4"),
        tickfont = list(size = 14)
      ),
      yaxis = list(
        title = list(text = "log fold change", font = list(size = 18, color = "black")),
        autorange = "reversed",
        tickvals = c(2, 1, 0, -1, -2),
        ticktext = c("2", "1", "0", "-1", "-2"),
        tickfont = list(size = 14)
      ),
      zaxis = list(
        title = list(text = "-log10(p-value)", font = list(size = 18, color = "black")),
        tickmode = "array",
        tickvals = c(2, 4, 6),
        ticktext = c("2", "4", "6"),
        tickfont = list(size = 14)
      ),
      yaxis = list(title = "log fold change"),
      zaxis = list(title = "-log10(p-value)"),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 0.9, z = 0.85),
      camera = list(eye = list(x = -0.75, y = -1.8, z = 0.7)
    )))


# === Save interactive plot ===
htmlwidgets::saveWidget(as_widget(phos3D_J),
                        "/Users/majaewen/Documents/Uni/6_Semester/Daten/Volcanos/PBMCs_PRO1.html")