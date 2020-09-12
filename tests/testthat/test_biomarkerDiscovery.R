context("biomarkerDiscovery")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

aucMVC_output <- aucMCV(cachexiaData, ref_level = "control")

test_that("aucMCV", {
    expect_equal(length(aucMVC_output), 12)
    
    expect_equal(as.character(aucMVC_output$ranking),
                 c("7.02077796211023", "5.47629546741624", "4.75284052850626", 
                   "4.55629155593845", "4.32815790692361", "3.49180858294745",
                   "3.49110785001989", "3.13744741886515"))
    
    expect_equal(names(aucMVC_output$ranking),
                 c("X3.Hydroxyisovalerate", "X3.Hydroxybutyrate", "X1.Methylnicotinamide",      
                   "X2.Hydroxyisobutyrate", "X3.Aminoisobutyrate", "X2.Oxoglutarate",            
                   "X2.Aminobutyrate", "X1.6.Anhydro.beta.D.glucose"))
})
