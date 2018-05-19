citHeader("To cite radjust in publications, please use:")

# citation(auto = meta)

year <- sub("-.*", "", meta$Date)
note <- sprintf("R package version %s", meta$Version)

bibentry(bibtype = "Manual",
         title = paste0("\"radjust: ",meta$Title,"\""),
         author =c(
           person("Shay", "Yaacoby", email = "shay66@gmail.com", role = "aut"),
           person("Marina", "Bogomolov", email = "marinabo@technion.ac.il", role = "aut"),
           person("Ruth", "Heller", email = "ruheller@gmail.com", role =  c("aut", "cre"))),
         year = year,
         note = note,
         url = "https://github.com/shay-y/radjust")


bibentry(bibtype = "article", # header = "to cite radjust_sym(), add:",
         author = c(person("Marina","Bogomolov"), person("Ruth","Heller")),
         title = "Assessing replicability of findings across two studies of multiple features",
         journal="arXiv preprint arXiv:1504.00534",
         year="2015"
)

bibentry(bibtype = "article", # header = "to cite radjust_pf(), add:",
         author = c(person("Ruth","Heller"), person("Marina","Bogomolov"), person("Yoav","Benjamini")),
         title = "Deciding whether follow-up studies have replicated findings in a preliminary large-scale omics study",
         journal = "Proceedings of the National Academy of Sciences",
         volume = "111",
         number = "46",
         pages = "16262--16267",
         year = "2014",
         publisher="National Acad Sciences"
)



