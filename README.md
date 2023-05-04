# Scripts and analyses for *R. hastatulus* XYY male genome paper

Initial proposed formatting for plots

```r
pubTheme <-
  theme(title = element_text(size=10), #usually 14
        text = element_text(size=8),
        plot.background = element_rect(fill="lightgrey"),
        strip.background = element_rect(linetype=0,linewidth=8,
                                        fill="grey"), #facet boxes
        strip.text = element_text(size=10), #usually 11
        legend.background = element_rect(fill="grey"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        panel.background = element_rect(fill="#FFFFFF"), #plot body
        panel.grid = element_blank() #plot body pt2
  )
pubColours <- c("#2d9da6","#99b700","#8a53b6","#bf4824","#152944")
```
