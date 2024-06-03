# Functions for uniform plotting and naming (implemented by Lucia)

long_names_pfts = function(x) {
  x = gsub("ibs", "Pioneering broadleaf", x)
  x = gsub("tebs", "Temperate broadleaf", x)
  x = gsub("bne", "Needleleaf evergreen", x)
  x = gsub("bine", "Shade-intolerant\nneedleleaf evergreen", x)
  x = gsub("bns", "Needleleaf summergreen", x)
  x = gsub("tene", "Temperate needleleaf", x)
  x = gsub("tundra", "Tundra", x)
  x = gsub("soil", "Bare soil", x)
  x = gsub("mixed forest", "Mixed forest", x)
  x = gsub("otherc", "Conifers (other)", x)
  x = gsub("regeneration failure", "Regeneration failure", x)
  return(x)
}


long_names_scenarios = function(x) {
  x = gsub(" control", " Historic", x)
  x = gsub("picontrol", "Control", x)
  x = gsub("ssp126", "SSP1-RCP2.6", x)
  x = gsub("ssp370", "SSP3-RCP7.0", x)
  x = gsub("ssp585", "SSP5-RCP8.5", x)
  
  return(x)
}


add_common_layout = function() {
  theme_classic() %+replace%
    theme(legend.background = element_rect(fill='transparent', color = NA), # make plot background transparent, especially helpful for presentations
          legend.box.background = element_rect(fill='transparent', color = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),  
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "#e5e5e5ff")) #zero margins to make paneling sub panels more contolled
}
