# these are some commands that tell this demo where to find the demo data
otu.table.filename <- system.file('demo', 'otus.dat', package='texmexseq')

# first, load up an OTU table of raw counts
# "row.names=1" makes the first column (the OTU IDs) become the row names
otu.table <- read.table(otu.table.filename, header=T, row.names=1)

# you can evaluate how good the texmex fit will be for any particular sample by making
# a PP plot. (it's a ggplot object, so you can add layers and manipulate the
# labels as you like.)
ppplot(otu.table$inoculum1.control.before)

# transform all the counts to F values. (you could transform to z values using
# the command z.transform.table instead.)
f.table <- f.transform.table(otu.table)

# now pull out a "quad", four before/after control/treatment samples.
# this command pulls out the four columns from the f.table, pulls out
# the OTU IDs (which were previously just rownames), and computes the
# change in F in the control and treatment
quad <- quad.table(f.table, 'inoculum1.control.before', 'inoculum1.control.after',
                            'inoculum1.treatment.before', 'inoculum1.treatment.after')

# a plot of this quad will compare the changes in F against one another
# like ppplot, quad.plot produces a ggplot object.
p <- quad.plot(quad)
p

# the plot might help you decide what your OTUs of interest will be.
# the filter command from dplyr can really help here.
# here, i'm picking all those OTUs that are in the top-left half-diamond
# of the treatment deltaF vs. control deltaF plot
interesting <- filter(quad, d.treatment > d.control + 1.0)

# i can add those in a separate color onto the plot
p + geom_points(data=interesting, col='red')

# and i can ask which OTUs are those interesting ones (with the ones that have 
# the biggest increase in F at the top). i use "head" to only look at the top
# 10 OTUs
arrange(interesting, desc(d.treatment)) %>% head
