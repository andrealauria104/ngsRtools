# Themes definitions ----
# Theme 1 --
my_theme <- theme(legend.position = "bottom"
                  , line = element_line(size=0.5)
                  , legend.title   = element_text(color="black",size=8)#rel(1.5)
                  , legend.text    = element_text(color="black",size=8)
                  , axis.title     = element_text(color="black",size=8)
                  , axis.text      = element_text(color="black",size=8)
                  , axis.line      = element_line(size=0.25)
                  , panel.grid.major = element_line(size=0.25)
                  , panel.grid.minor = element_blank()
                  , strip.text       = element_text(size=8)
                  , strip.background = element_rect(fill = NA)
                  , plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))
# Theme 2 --
my_theme_2 <- theme(legend.position = "bottom"
                  , line = element_line(size=0.5)
                  , legend.title   = element_text(color="black",size=8)#rel(1.5)
                  , legend.text    = element_text(color="black",size=8)
                  , axis.title     = element_text(color="black",size=8)
                  , axis.text      = element_text(color="black",size=8)
                  , axis.line      = element_line(size=0.25)
                  , panel.grid.major = element_line(size=0.25)
                  , panel.grid.minor = element_blank()
                  , strip.text       = element_text(size=8)
                  , strip.background = element_rect(fill = NA)) +
  theme(strip.background = element_blank()
        , strip.text = element_text(size = 8)
        , text = element_text(size = 8)
        , plot.title = element_text(size = 8, hjust = 0.5, face = "bold")
        , plot.background = element_rect(size = 0.25)
        , panel.background = element_rect(size = 0.25)
        , panel.border = element_rect(size=0.25)
        , panel.grid = element_line(size = 0.25)
        , axis.ticks = element_line(size = 0.25))

# Theme blank --
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank()
    , axis.title.y = element_blank()
    , panel.border = element_blank()
    , panel.grid=element_blank()
    , axis.ticks = element_blank()
    , plot.title=element_text(size=8, face="bold")
    , line = element_line(size=0.5)
    , axis.title     = element_text(color="black",size=8)
    , axis.text      = element_text(color="black",size=8)
    , legend.title   = element_text(color="black",size=8)#rel(1.5)
    , legend.text    = element_text(color="black",size=8)
  )
