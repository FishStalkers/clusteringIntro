pacman::p_load(bnlearn,Rgraphviz, gRain)

asia <- bnlearn::asia
colnames(asia) <- c("VisittoAsia","Smoking","Tuberculosis","LungCancer","Bronchitis", "TuberculosisCersusLungCancer_Bronchitis",
                      "ChestXray","Dyspnoea")

bn_df <- data.frame(asia)
res <- hc(bn_df)
plot(res)
res = model2network("[VisittoAsia][Smoking][Tuberculosis|VisittoAsia][LungCancer|Smoking][Bronchitis|Smoking][Dyspnoea|Bronchitis:TuberculosisCersusLungCancer_Bronchitis][TuberculosisCersusLungCancer_Bronchitis|Tuberculosis:LungCancer][ChestXray|TuberculosisCersusLungCancer_Bronchitis]")
plot(res)
fittedbn <- bn.fit(res, data = bn_df)
graphviz.chart(
  fittedbn,
  type = "barprob" ,
  layout = "fdp",
  scale = c(0.75, 0.9),
  grid = TRUE,
  bar.col = "red",
  strip.bg = "lightskyblue"
)
fittedbn <- bn.fit(res, data = bn_df)
print(fittedbn$Tuberculosis)
print(fittedbn$TuberculosisCersusLungCancer_Bronchitis)
cpquery(fittedbn, event = (Smoking=="yes"), evidence = (Bronchitis=="no"))
cpquery(fittedbn, event = (Smoking=="yes"), evidence = (Bronchitis=="yes"))
cpquery(fittedbn, event = (Tuberculosis=="yes"), evidence = (Dyspnoea=="yes"))
cpquery(fittedbn, event = (Tuberculosis=="no"), evidence = (Dyspnoea=="yes"))


