PrintRunNumbers()
{
  string runType;
  int	 runNumber;
  int	 goodRunsCounter = 0;

  ofstream outfile;
  outfile.open("Good_LED_runNumbers.txt", ios::app);

  for(int i = 80; i < 122; i++)
  {
//    if(i == 91 || i == 93 || i == 101 || i == 107 || i == 121)
    if((i >= 82 && i <= 87) || (i >= 90 && i <= 93) || (i >= 96 && i <= 101)
	|| i == 107 || (i >= 110 && i <= 113) || i == 115 || i == 117 || i == 121)
    {
      continue;
    }

    TString fileName = Form("OctetLists/octet_list_%i.dat", i);

    string buf1;
    ifstream infile1;
    cout << "The file being opened is: " << fileName << endl;
    infile1.open(fileName);

    //a check to make sure the file is open
    if(!infile1.is_open())
      cout << "Problem opening " << fileName << endl;

    while(true)
    {
      getline(infile1, buf1);
      istringstream bufstream1(buf1);

      if(!infile1.eof())
      {
        bufstream1 >> runType >> runNumber;
      }

      if(infile1.eof())
        break;

      if(runNumber == 21792 || runNumber == 21903 || runNumber == 21932 || runNumber == 21933
	|| runNumber == 21777 || runNumber == 21891 || runNumber == 21948 || runNumber == 21911
	|| runNumber == 22261	// the previous runs are all "cutruns" according to Simon
	|| (runNumber >= 21827 && runNumber <= 21842) || (runNumber >= 21865 && runNumber <= 21877)
	|| (runNumber >= 21993 && runNumber <= 22003) || (runNumber >= 22084 && runNumber <= 22093)
	|| (runNumber >= 22123 && runNumber <= 22163) || (runNumber >= 22311 && runNumber <= 22339)
	|| (runNumber >= 22344 && runNumber <= 22350) || (runNumber >= 22358 && runNumber <= 22389)
	|| (runNumber >= 22412 && runNumber <= 22437) || (runNumber >= 22441 && runNumber <= 22480)
	|| (runNumber >= 22483 && runNumber <= 22511) || (runNumber >= 22733 && runNumber <= 22737)
	|| (runNumber >= 22746 && runNumber <= 22753) || (runNumber >= 22793 && runNumber <= 22803)
	|| (runNumber >= 22829 && runNumber <= 22834) || (runNumber >= 22873 && runNumber <= 22896)
	|| (runNumber >= 22947 && runNumber <= 22953) || (runNumber >= 22955 && runNumber <= 22980)
	|| (runNumber >= 23010 && runNumber <= 23017) || (runNumber >= 23081 && runNumber <= 23084)
	// the above cut runs are all called bimodal ranges
	|| (runNumber >= 22453 && runNumber <= 22463) || runNumber == 22222 || runNumber == 22301
	|| runNumber == 22302 || runNumber == 22448 || runNumber == 22450
	|| runNumber == 22778 // the previous cut runs are all called ledcullrange
        )
     {
       continue;
     }


      outfile << i << "\t"
	      << runType << "\t"
	      << runNumber << "\n";

      goodRunsCounter++;
    }


    if(goodRunsCounter >= 24)
    {
      cout << "Octet " << i << " has " << goodRunsCounter << " good runs." << endl;
    }

    goodRunsCounter = 0;

  }

  outfile.close();

  cout << "Finished printing run numbers to Good_LED_runNumbers.txt" << endl;

}
