using ExcelReaders

function generatedCSVfiles(file::String)
   println("Processing parameter B... 1/4")
   B = readxlsheet(file, "B" )
   writecsv("B.csv", B)

   println("Processing parameter NV... 2/4")
   NV = readxlsheet(file, "NV" )
   writecsv("NV.csv", NV)

   println("Reading parameter ST... 3/4")
   ST = readxlsheet(file, "ST" )
   writecsv("ST.csv", ST)

   println("Processing parameter NV... 4/4. This one might takes a while...")
   TT = readxlsheet(file, "TT" )
   writecsv("TT.csv", TT)
end


function readData(totalN, totalJ, totalV, totalS)

   println("Reading parameter B... 1/4")
   B = readcsv("B.csv" )
   B = B[:, 1:totalS]

   #EC = zeros(totalJ)
   #for j = 1:totalJ EC[j] = rand(50:100) end
   EC = 0.+[5 10 15 20]

   M = 500

   println("Reading parameter NV... 2/4")
   NV = readcsv("NV.csv")
   NV = NV[:, 1:totalS]

   println("Reading parameter ST... 3/4")
   ST = readcsv("ST.csv")
   ST = ST[:, 1:totalS]

   P = 1/totalS.*ones(totalS)

   println("Reading parameter TT... 4/4")
   TTaux = readcsv("TT.csv")
   TTaux= reshape(TTaux[:, 1:totalS], (totalN,totalN,1,totalS))

   TT = Array(Float64, (totalN, totalN, totalV, totalS))

   for v in 1:totalV
      TT[:, :, v, :] = 1.*TTaux[:, :, 1, :]
   end

   println("Done.")

   return B, EC, M, NV, ST, P, TT
end
