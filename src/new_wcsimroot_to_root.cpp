#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom1.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"

using namespace std;




void get_signal_from_WCSim(string infilename, string outfilename, const map<string, string> config, float event_info[2], int max_events) {
    
    
    
    // Extracted parameters from config 
    //const int wcsim_version = std::stoi(config.at("WCSIM_VERSION")); // Should be 0 or 1 for now - Erwan 11/12/2024

    const int maxHitsSig = std::stoi(config.at("MAX_HITS_SIG"));
    const short int verbose = std::stoi(config.at("VERBOSE"));
    Int_t eventType = (Int_t) std::stoi(config.at("EVENT_TYPE"));
    
    //const int maxHitsBg = std::stoi(config["MAX_HITS_BG"]);
    //const double minDurationBg = std::stod(config["MIN_DURATION_BG"]);



    /*************************************************************************/
    // ####  Get global defintion of the data from the input files
    /*************************************************************************/

    // Open the file containing the tree, and its geometry
    TFile *file    = TFile::Open(infilename.data(), "READ");

    // Get the wcsimgeoT TTree object from the input file (geo for geometry)
    TTree *geotree = (TTree*) file->Get("wcsimGeoT");

    WCSimRootGeom *geo = new WCSimRootGeom();        // Special trees of wcsim requier extra data manager. Classe pour GÉRER le TTree geotree (compatibilité root & WCSIM)
    geotree->SetBranchAddress("wcsimrootgeom", &geo);


    // Get the wcsimT TTree object from the input file
    // Use of the WCSimRoot special classES to deal with root & WCSim
    TTree *wcsimT= (TTree*)file -> Get("wcsimT"); 	
    WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();  // Special trees of wcsim requier extra data manager Create a WCSimRootEvent to store data from the tree in
    wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);     // Force deletion to prevent memory leak 


    TBranch *branch = wcsimT->GetBranch("wcsimrootevent");
    branch->SetAddress(&wcsimrootsuperevent);


    // Test if the inputfile isn't empty
    if( geotree->GetEntries() == 0 ){ 
        std::cout << "ERROR, there is 0 event in the Geometry TTree" << endl; 
        exit(2); 
    }


    /*************************************************************************/
    // ####  Prepare data for output file
    /*************************************************************************/

    // --- Instantiate the variable destinated to the tree branches
    // Event level
    Int_t n_hits;                                           

    Float_t energy;                                                        
    Float_t vertex[3], particleDir[3], particleStart[3], particleStop[3];   
    Float_t vertex_time, time_trigger;

    // PMT level
    Int_t tubeIds[maxHitsSig];                           
    Float_t charge[maxHitsSig];
    Float_t time[maxHitsSig];       
    Float_t hitx[maxHitsSig], hity[maxHitsSig], hitz[maxHitsSig];

    // Creates the new file and its tree
    TFile *outfile = new TFile(outfilename.data(), "RECREATE");
    TTree *outtree = new TTree("pure_root_tree", "events_data");

    // --- Create the output tree branches
    // Event level
    outtree->Branch("eventType", &eventType, "eventType/I");

    outtree->Branch("energy", &energy, "energy/F");
    outtree->Branch("vertex", &vertex, "vertex[3]/F");
    outtree->Branch("vertex_time", &vertex_time, "vertex_time/F");

    outtree->Branch("particleDir", &particleDir, "particleDir[3]/F");
    outtree->Branch("particleStart", &particleStart, "particleStart[3]/F");
    outtree->Branch("particleStop", &particleStop, "particleStop[3]/F");

    outtree->Branch("n_hits", &n_hits, "n_hits/I");
    outtree->Branch("time_trigger", &time_trigger, "time_trigger/F");


    // PMT level
    outtree->Branch("tubeIds", &tubeIds, "tubeIds[n_hits]/I"); 
    outtree->Branch("hitx", &hitx, "hitx[n_hits]/F");
    outtree->Branch("hity", &hity, "hity[n_hits]/F");
    outtree->Branch("hitz", &hitz, "hitz[n_hits]/F");
    outtree->Branch("charge", &charge, "charge[n_hits]/F");
    outtree->Branch("time", &time, "time[n_hits]/F");


    /*************************************************************************/
    // Get high-level characteristics of the data stored in the geoT
    /*************************************************************************/
    float cylinderRadius, cylinderHeight;
    float PMTradius[2];
    int nevents = 0;

    geotree->GetEntry(0);
    nevents = wcsimT->GetEntries();	
    if (verbose >= 1) {

        cout << "\n\n Parsing geometry data.. \n\n"  << endl;


        cylinderRadius =  geo -> GetWCCylRadius();
        cylinderHeight = geo -> GetWCCylLength();
        cout << "Cylinder radius : " << cylinderRadius << endl;
        cout << "Cylinder height : " << cylinderHeight << endl;

        PMTradius[0]=geo->GetWCPMTRadius();
        PMTradius[1]=geo->GetWCPMTRadius(true);
        cout << "\nNumber of PMTs of 1st type = " << geo->GetWCNumPMT() << ", radius = " << PMTradius[0] << endl;
        cout << "Number of PMTs of 2nd type = " << geo->GetWCNumPMT(true) << ", radius = " << PMTradius[1] <<endl;

        cout << "\n -- Looking for data in the file : " << infilename << endl;
        cout << " -- The file contains " << nevents << " events  -- \n" << endl;

        if ( verbose >= 2) { cout << " -- Event type : " << eventType << "\n" << endl;}

        cout << "Finished parsing geoTree. \nStarting to process events.."  << endl;
    }


    // Determine the range of events to process
    int end_event = (max_events > 0 && max_events < nevents) ? max_events : nevents;

    if (verbose >= 1) {
        std::cout << "Processing events 0 to " << end_event - 1 << " out of " << nevents << " total events." << endl;
    }
    
    /*************************************************************************/
    // Loop over the EVENTS
    /*************************************************************************/
    // Variables for event management

    // Pointors for wcsimevent
    WCSimRootTrigger *wcsimrootevent = nullptr;

    // Pointors for GetNTrack
    TObject *element =  nullptr;
    WCSimRootTrack *wcsimroottrack = nullptr;

    // Pointors for GetNcherenkovdigihits
    TObject *Hit = nullptr;
    WCSimRootCherenkovDigiHit *cDigiHit = nullptr;
    WCSimRootPMT hitpmt; 

    std::vector<double> triggerInfo; // vector<double> class type mandatory by wcsim
    double triggerShift, triggerTime;


    int ntrack = 0;
    int nevents_check = 0;

    for ( int i = 0; i < end_event; i++ ) {
         
        wcsimT->GetEvent(i);                                     //  Load the i^th event into wcsimrootsuperevent
        wcsimrootevent = wcsimrootsuperevent -> GetTrigger(0);   // load the trigger of the i_th event 

        triggerInfo = wcsimrootevent->GetTriggerInfo();
        if ( triggerInfo.size() >= 3 ) {
            triggerShift = triggerInfo[1];
            triggerTime = triggerInfo[2];

            if ( triggerShift != 950.0 ) {
                std::cout << "\n\nUnexpected value on triggerShift \nEvent num" << i << endl;
                std::cout << "triggerInfo size " << triggerInfo.size()  << endl;
                std::cout << "N_hits if this event : " << wcsimrootevent->GetNcherenkovdigihits() << endl;
                std::cout << "triggerShift " << triggerShift  << endl;
                std::cout << "triggerTime " << triggerTime  << endl;
                break;
            }
        } else { triggerTime = -1.; } // Sentinel value for events with no hits }

        time_trigger = triggerTime;


        // Get the vertex of the mother particle
        for (int k = 0; k < 3; k++) { 
            vertex[k] = (Float_t) wcsimrootevent ->GetVtx_old(k);
        }

        // Ntrack deals with the data about the mother particule. 
        ntrack = wcsimrootevent->GetNtrack();
        for ( int k=0 ; k < ntrack ; k++ ) {

            element = (wcsimrootevent->GetTracks())->At(k);
            wcsimroottrack = dynamic_cast<WCSimRootTrack*> (element);

            if (wcsimroottrack->GetIpnu() != 0 && wcsimroottrack->GetFlag() == -1) {
                
                energy = wcsimroottrack->GetE();
                vertex_time = wcsimroottrack->GetTime();
                for (int k = 0; k < 3; k++) {
                    particleDir[k]   = wcsimroottrack->GetDir(k);
                    particleStart[k] = wcsimroottrack->GetStart(k);
                    particleStop[k]  = wcsimroottrack->GetStop(k);
                }
                
                break; 
            // C'est quoi Ipnu, c'est quoi Getflag, c'est quoi le but de cette ligne (et de cette boucle for)
            // Hypothèse : y'a qu'un seul k qui vérifie cette condition, si on l'a trouvé on peut sortir de la boucle for
            // Mais pourquoi la particule mère ne se trouve pas au tout début de la track ? Donc a k = 0 ou k = n_track - 1
            }
        }

        // Loop over the HITS of one event
        n_hits = wcsimrootevent->GetNcherenkovdigihits();

        for(int j=0 ; j<n_hits ; j++) {

            Hit = ( wcsimrootevent->GetCherenkovDigiHits() ) -> At(j);
            cDigiHit = dynamic_cast<WCSimRootCherenkovDigiHit*> (Hit);
            hitpmt = geo->GetPMT(cDigiHit -> GetTubeId() - 1, false); // Pourquoi - 1 ? Est-ce que one_indexe=False dans watchaml alors ?

            charge[j]  = cDigiHit->GetQ();
            time[j]    = cDigiHit->GetT();
            hitx[j]    = hitpmt.GetPosition(0);
            hity[j]    = hitpmt.GetPosition(1);
            hitz[j]    = hitpmt.GetPosition(2);
            tubeIds[j] = cDigiHit -> GetTubeId();

        }

        nevents_check++;
        outtree -> Fill();

        // Loop print
        if ( (i % 1000 == 0) && (verbose >= 1) ) {
            cout << "Event #" << i << ", found " << n_hits << " Cherenkov hits" << endl;

            if ( verbose >= 2) {
            cout << "Event #" << i << ", ID=" << wcsimroottrack->GetIpnu() << "  Pnu=" << energy << endl;
            }
        }
    }

    outtree->Write();
    outfile->Close();

    if ( verbose >=1 ) {
        cout << "\n -- End of the processing -- \n" << endl;

        cout << "There is " << nevents       << " events in the input file." << endl;
        cout << "There is " << nevents_check << " events in the output file." << endl;
        cout << "\n -- Saving the data at : " << outfilename << endl;
    }

    // Update event_info with results
    event_info[0] = -1.; // No idea of the purpose of this
    event_info[1] = static_cast<float>(nevents);

}


int main(int argc, char* argv[]) {

    if (argc < 2 || argc > 3) {
        cerr << "\n --- Error --- \n Usage: " << argv[0] << " path_to_config.txt [max_events] \n" << endl;
        return 1;
    }

    std::ifstream configFile(argv[1]);
    std::map<std::string, std::string> configValues;

    std::string line;
    while (std::getline(configFile, line)) {
        // Check if the line starts with '#' and skip it if it does
        if (!line.empty() && line[0] == '#') { continue; }

        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                configValues[key] = value;
            }
        }
    }

    // Parameters for get_signal_from_WCSim
    const string input_file_path  = configValues["WCSIMROOT_FILE_PATH"];
    const string output_file_path = configValues["RESULT_FILE_PATH"];

    // Optional max events parameter
    int max_events = -1; // Default: process all events
    if (argc == 3) {
        max_events = std::stoi(argv[2]);
    }

    // Declare the float array for results
    float event_info[2] = {0.0f, 0.0f};

    // Call get_signal_from_WCSim
    get_signal_from_WCSim(input_file_path, output_file_path, configValues, event_info, max_events);

    // Print results
    //cout << "Total Event Duration: " << event_info[0] << endl;
    cout << "Number of Events Processed: " << event_info[1] << endl;

    return 0;
}
