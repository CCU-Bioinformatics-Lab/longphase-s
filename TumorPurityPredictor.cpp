#include "TumorPurityPredictor.h"

TumorPurityPredictor::TumorPurityPredictor(
    const HaplotagParameters& params,
    const std::vector<std::string>& chrVec,
    std::map<std::string, std::map<int, PosBase>>& chrPosNorBase,
    std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo
) : params(params), chrVec(chrVec), chrPosNorBase(chrPosNorBase), chrPosSomaticInfo(chrPosSomaticInfo), initial_data_size(0){}

TumorPurityPredictor::~TumorPurityPredictor(){

}

double TumorPurityPredictor::predictTumorPurity(){
   std::cerr << "predicting tumor purity... ";
    std::time_t begin = time(NULL);

    // median, iqr
    std::vector<PurityData> purityFeatureValueVec;
    double purity = 0.0;

    try{
        // build the feature value vector with first stage filter
        buildPurityFeatureValueVec(purityFeatureValueVec);

        std::cerr << "initial data size: " << initial_data_size << std::endl;
        std::cerr << "\n==========first filter==========" << std::endl;
        std::cerr << "[INFO] Data count: " << purityFeatureValueVec.size() << std::endl;
        std::cerr << "[INFO] consistencyRatioInNorBam == 0.0: " << filterCounts.consistencyRatioInNorBam << std::endl;
        std::cerr << "[INFO] consistencyRatio == 0.0: " << filterCounts.consistencyRatio << std::endl;
        std::cerr << "[INFO] consistencyRatioInNorBam <= 0.7: " << filterCounts.consistencyRatioInNorBamMaxThr << std::endl;
        std::cerr << "[INFO] readHpCountInNorBam <= 5: " << filterCounts.readHpCountInNorBam << std::endl;
        std::cerr << "[INFO] percentageOfGermlineHpInNorBam: >= 0.7 " << filterCounts.percentageOfGermlineHp << std::endl;

        // find the threshold of germlineReadHpCount for peak valley filter
        int germlineReadHpCountThreshold = findPeakValleythreshold(params, purityFeatureValueVec);

        // peak valley filter
        peakValleyFilter(purityFeatureValueVec, germlineReadHpCountThreshold);

        std::cerr << "\n==========second filter==========" << std::endl;
        std::cerr << "[INFO] germlineReadHpCountThreshold: " << germlineReadHpCountThreshold << std::endl;
        std::cerr << "[INFO] filter peakValley count: " << filterCounts.peakValley << std::endl;
        std::cerr << "[INFO] Data count: " << purityFeatureValueVec.size() << std::endl;

        BoxPlotValue plotValue = statisticPurityData(purityFeatureValueVec);

        //remove outliers for reduce noise 
        size_t iteration_times = 1;
        for(size_t i = 0; i < iteration_times; i++){
            std::printf("[INFO]==========iteration %ld==========\n", i+1);
            // remove outliers for reduce noise
            removeOutliers(purityFeatureValueVec, plotValue);
            // statistic the data
            plotValue = statisticPurityData(purityFeatureValueVec);
            std::printf("[INFO] wisker : lower = %f, upper = %f\n", plotValue.lowerWhisker, plotValue.upperWhisker);
            std::printf("[INFO] after remove outliers: %ld, data size: %ld\n", filterCounts.outliers, plotValue.data_size);
        }


        double median = plotValue.median;
        double iqr = plotValue.iqr;

        // purity prediction model
        purity = -3.3454 * median + 14.7747 * iqr + 4.0344 * median * median + -13.7777 * median * iqr + -5.2434 * iqr * iqr + 0.3058;

        if(purity > 1.0){
            purity = 1.0;
        }else if(purity < 0.0){
            throw std::runtime_error("The value of purity exceeds the model's prediction range: " + std::to_string(purity));
        }

        std::cerr<< difftime(time(NULL), begin) << "s\n";
        std::cerr << "[INFO] ===== tumor purity ===== : " << purity << std::endl;
        std::cerr << "[INFO] median: " << plotValue.median << std::endl;
        std::cerr << "[INFO] iqr: " << plotValue.iqr << std::endl;
        std::cerr << "[INFO] q1: " << plotValue.q1 << std::endl;
        std::cerr << "[INFO] q3: " << plotValue.q3 << std::endl;
        std::cerr << "[INFO] lowerWhisker: " << plotValue.lowerWhisker << std::endl;
        std::cerr << "[INFO] upperWhisker: " << plotValue.upperWhisker << std::endl;
        std::cerr << "[INFO] outliers: " << plotValue.outliers << std::endl;
        
        // write purity log
        writePurityLog(params, purity, plotValue, iteration_times, germlineReadHpCountThreshold);
    }catch(const std::exception& e){
        std::cerr << "[ERROR] " << e.what() << std::endl;
        std::cerr << "[ERROR] Failed to predict tumor purity, set purity to 0.0" << std::endl;
        purity = 0.0;
    }

    return purity;
}

void TumorPurityPredictor::buildPurityFeatureValueVec(std::vector<PurityData> &purityFeatureValueVec){
    for(auto chr : chrVec){
        std::map<int, HP3_Info>::iterator somaticPosIter = chrPosSomaticInfo[chr].begin();
        while(somaticPosIter != chrPosSomaticInfo[chr].end()){
            // int H1readCount = (*somaticPosIter).second.base.ReadHpCount[ReadHP::H1];
            // int H2readCount = (*somaticPosIter).second.base.ReadHpCount[ReadHP::H2];
            int pos = (*somaticPosIter).first;
            initial_data_size++;

            // the ratio of the read count of H1 and H2
            double germlineReadHpImbalanceRatio = (*somaticPosIter).second.base.germlineHaplotypeImbalanceRatio;

            //read hp count in the normal bam
            int H1readCountInNorBam = chrPosNorBase[chr][pos].ReadHpCount[ReadHP::H1];
            int H2readCountInNorBam = chrPosNorBase[chr][pos].ReadHpCount[ReadHP::H2];
            int germlineReadHpCountInNorBam = H1readCountInNorBam + H2readCountInNorBam;

            double germlineReadHpImbalanceRatioInNorBam = chrPosNorBase[chr][pos].germlineHaplotypeImbalanceRatio;

            double percentageOfGermlineHpInNorBam = chrPosNorBase[chr][pos].percentageOfGermlineHp;

            bool includeInStatistics = true;

            // statistic the filter out data
            if(germlineReadHpImbalanceRatioInNorBam == GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MIN_THR){
                includeInStatistics = false;
                filterCounts.consistencyRatioInNorBam++;
            }else if(germlineReadHpImbalanceRatio == GERMLINE_HP_IMBALANCE_RATIO_MIN_THR){
                includeInStatistics = false;
                filterCounts.consistencyRatio++;
            }else if(germlineReadHpImbalanceRatioInNorBam >= GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MAX_THR){
                includeInStatistics = false;
                filterCounts.consistencyRatioInNorBamMaxThr++;
            }else if(germlineReadHpCountInNorBam <= GERMLINE_HP_READ_COUNT_IN_NOR_BAM_MIN_THR){
                includeInStatistics = false;
                filterCounts.readHpCountInNorBam++;
            }else if(percentageOfGermlineHpInNorBam <= GERMLINE_HP_PERCENTAGE_IN_NOR_BAM_MAX_THR){
                includeInStatistics = false;
                filterCounts.percentageOfGermlineHp++;
            }else if(includeInStatistics){
                purityFeatureValueVec.emplace_back(
                    PurityData{
                        .chr = chr,
                        .pos = (*somaticPosIter).first,
                        .germlineReadHpConsistencyRatio = germlineReadHpImbalanceRatio,
                        .germlineReadHpCountInNorBam = germlineReadHpCountInNorBam
                    }
                );
                // (*somaticPosIter).second.statisticPurity = true;
                chrPosSomaticFlag[chr][(*somaticPosIter).first].statisticPurity = true;
            }
            somaticPosIter++;
        }
    }

    if(purityFeatureValueVec.empty()) {
        throw std::runtime_error("Failed to build purity feature vector: empty vector");
    }
}


int TumorPurityPredictor::findPeakValleythreshold(const HaplotagParameters& params, const std::vector<PurityData> &purityFeatureValueVec){
    
    std::cerr << "[INFO] peak valley filter ..." << std::endl;
    int threshold = 0;

    try{
        Histogram histogram;
        // build histogram
        histogram.buildHistogram(purityFeatureValueVec);
        // calculate the percentage of data and resize the histogram
        histogram.calculateStatistics();

        double sigma = 0.5;
        std::cerr << "[INFO] apply gaussian filter with sigma: " << sigma << std::endl;
        // apply gaussian filter
        Histogram smoothedHistogram = histogram.getSmoothedHistogram(sigma);

        PeakProcessor peakSet;

        static constexpr double MIN_PEAK_RATIO = 0.05;
        static constexpr int MIN_DISTANCE = 2;

        double max_height = smoothedHistogram.getMaxHeight();
        size_t total_snp_count = smoothedHistogram.getTotalSnpCount();
        std::pair<size_t, size_t> data_range = smoothedHistogram.getDataRange();

        double peak_threshold = std::max(
            static_cast<size_t>((double)max_height * MIN_PEAK_RATIO), 
            static_cast<size_t>(1)
        );

        //find the peaks with the threshold
        peakSet.findPeaks(smoothedHistogram.getHistogram(), peak_threshold);
        //remove close peak with min_distance
        peakSet.removeClosePeaks(MIN_DISTANCE);
        //determine the trend of the peak
        peakSet.determineTrends();
        //find the main peak
        peakSet.findMainPeakCandidates();
        //set the threshold by the lowest valley
        peakSet.SetThresholdByValley(smoothedHistogram.getHistogram(), smoothedHistogram.getMaxHeight());

        threshold = peakSet.getThreshold();

        // print the exec_log
        for(const auto& log : peakSet.exec_log){
            std::cout << log << std::endl;
        }

        //write histogram to file
        if(params.writeReadLog){
            peakSet.writePeakValleyLog(params, 
                                histogram.getHistogram(),  
                                smoothedHistogram.getHistogram(),  
                                total_snp_count, 
                                data_range, 
                                max_height, 
                                MIN_PEAK_RATIO, 
                                peak_threshold,
                                sigma);
        }

    }catch(const std::exception& e){
        std::cerr << "[ERROR] " << e.what() << std::endl;
        std::cerr << "[ERROR] Failed to find peak valley threshold, set threshold to 0" << std::endl;
        threshold = 0;
    }
    
    return threshold;
}

void TumorPurityPredictor::peakValleyFilter(std::vector<PurityData> &purityFeatureValueVec, int &germlineReadHpCountThreshold){
    std::vector<PurityData>::iterator vecFeatureIter = purityFeatureValueVec.begin();
    while(vecFeatureIter != purityFeatureValueVec.end()){
        if((*vecFeatureIter).germlineReadHpCountInNorBam < germlineReadHpCountThreshold){
            filterCounts.peakValley++;
            // flag the position as not used for purity prediction
            // chrPosSomaticInfo[(*vecFeatureIter).chr][(*vecFeatureIter).pos].statisticPurity = false;
            chrPosSomaticFlag[(*vecFeatureIter).chr][(*vecFeatureIter).pos].statisticPurity = false;
            vecFeatureIter = purityFeatureValueVec.erase(vecFeatureIter);
        }else{
            vecFeatureIter++;
        }
    }
}

void TumorPurityPredictor::removeOutliers(std::vector<PurityData> &purityFeatureValueVec, BoxPlotValue &plotValue){
     std::vector<PurityData>::iterator vecIter = purityFeatureValueVec.begin();

    while(vecIter != purityFeatureValueVec.end()){
        if((*vecIter).germlineReadHpConsistencyRatio < plotValue.lowerWhisker || (*vecIter).germlineReadHpConsistencyRatio > plotValue.upperWhisker){
            //remove outliers for reduce noise
            std::string chr = (*vecIter).chr;
            int pos = (*vecIter).pos;
            //flag the position as not used for purity prediction
            // chrPosSomaticInfo[chr][pos].statisticPurity = false;
            chrPosSomaticFlag[chr][pos].statisticPurity = false;
            vecIter = purityFeatureValueVec.erase(vecIter);
            filterCounts.outliers++;
        }else{
            vecIter++;
        }
    }
}
BoxPlotValue TumorPurityPredictor::statisticPurityData(std::vector<PurityData> &purityFeatureValueVec){
    BoxPlotValue plotValue;
    plotValue.data_size = purityFeatureValueVec.size();
    
    try{
        // Handle empty data case
        if (plotValue.data_size == 0) {
            throw std::runtime_error("the data size is 0");
        }

        // Sort the original data
        std::sort(purityFeatureValueVec.begin(), purityFeatureValueVec.end(), PurityData::compareByGermlineReadHpConsisRatio);

        // Define percentile calculation function using linear interpolation
        auto percentile = [&](double p) -> double {
            if (p < 0.0 || p > 1.0) {
                throw std::invalid_argument("Percentile must be between 0 and 1, p = " + std::to_string(p));
            }
            
            double pos = p * (plotValue.data_size - 1);
            size_t idx = static_cast<size_t>(pos);
            double frac = pos - idx;
            
            //check index 
            if (idx + 1 >= plotValue.data_size) {
                return purityFeatureValueVec[plotValue.data_size - 1].germlineReadHpConsistencyRatio;
            }

            // //element is the last data
            if (idx == plotValue.data_size - 1) {
                return purityFeatureValueVec[idx].germlineReadHpConsistencyRatio;
            }
            
            // return linear interplolation
            return purityFeatureValueVec[idx].germlineReadHpConsistencyRatio * (1.0 - frac) + 
                purityFeatureValueVec[idx + 1].germlineReadHpConsistencyRatio * frac;
        };

        // Calculate quartiles and median
        plotValue.q1 = percentile(0.25);      // First quartile
        plotValue.median = percentile(0.5);    // Median
        plotValue.q3 = percentile(0.75);      // Third quartile
        plotValue.iqr = plotValue.q3 - plotValue.q1;  // Interquartile range
        
        // Calculate whiskers (1.5 * IQR rule)
        plotValue.lowerWhisker = std::max(0.0, plotValue.q1 - 1.5 * plotValue.iqr);
        plotValue.upperWhisker = plotValue.q3 + 1.5 * plotValue.iqr;
        
        // Count outliers (points beyond whiskers)
        plotValue.outliers = 0;
        for (const auto& value : purityFeatureValueVec) {
            if (value.germlineReadHpConsistencyRatio < plotValue.lowerWhisker || value.germlineReadHpConsistencyRatio > plotValue.upperWhisker) {
                plotValue.outliers++;
            }
        }
            
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Invalid argument in statistics calculation: " + std::string(e.what()));
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to statistic purity data: " + std::string(e.what()));
    }

    return plotValue;
}

void TumorPurityPredictor::markStatisticFlag(std::map<std::string, std::map<int, HP3_Info>>& chrPosSomaticInfo){
    for(const auto& chr : chrVec) {
        auto& chrInfo = chrPosSomaticInfo[chr];
        for(const auto& pos : chrPosSomaticFlag[chr]) {
            auto it = chrInfo.find(pos.first);
            if (it != chrInfo.end()) {
                it->second.statisticPurity = pos.second.statisticPurity;
            } else {
                std::cerr << "[WARNING] Faill to flag statistic purity for the position : chr: " << chr << " pos: " << pos.first << std::endl;
            }
        }
    }
}

void TumorPurityPredictor::writePurityLog(const HaplotagParameters &params, double &purity, BoxPlotValue &plotValue, size_t &iteration_times, int &germlineReadHpCountThreshold){
    try{
        std::ofstream purityLog = std::ofstream(params.resultPrefix+"_purity.out");

        if(!purityLog.is_open()){
            throw std::runtime_error("Failed to open purity log file: " + params.resultPrefix+"_purity.out");
        }

        purityLog << "#Initial data size: " << initial_data_size << std::endl;
        purityLog << "#==========filter parameters==========" << std::endl;
        purityLog << "#GERMLINE_HP_CONSISTENCY_RATIO_MIN_THR: " << GERMLINE_HP_IMBALANCE_RATIO_MIN_THR << std::endl;
        purityLog << "#GERMLINE_HP_CONSISTENCY_RATIO_IN_NOR_BAM_MIN_THR: " << GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MIN_THR << std::endl;
        purityLog << "#GERMLINE_HP_CONSISTENCY_RATIO_IN_NOR_BAM_MAX_THR: " << GERMLINE_HP_IMBALANCE_RATIO_IN_NOR_BAM_MAX_THR << std::endl;
        purityLog << "#GERMLINE_HP_PERCENTAGE_IN_NOR_BAM_MAX_THR: " << GERMLINE_HP_PERCENTAGE_IN_NOR_BAM_MAX_THR << std::endl;
        purityLog << "#GERMLINE_HP_READ_COUNT_IN_NOR_BAM_MIN_THR: " << GERMLINE_HP_READ_COUNT_IN_NOR_BAM_MIN_THR << std::endl;
        purityLog << "#GERMLINE_HP_READ_COUNT_IN_NOR_BAM_DYNAMIC_THR: " << germlineReadHpCountThreshold << std::endl;

        purityLog << "#==========initial filter out data count==========" << std::endl;
        purityLog << "#consistencyRatioInNorBam: " << filterCounts.consistencyRatioInNorBam << std::endl;
        purityLog << "#consistencyRatio: " << filterCounts.consistencyRatio << std::endl;
        purityLog << "#consistencyRatioInNorBam_over_thr: " << filterCounts.consistencyRatioInNorBamMaxThr << std::endl;
        purityLog << "#readHpCountInNorBam: " << filterCounts.readHpCountInNorBam << std::endl;
        purityLog << "#percentageOfGermlineHpInNorBam: " << filterCounts.percentageOfGermlineHp << std::endl;

        purityLog << "#==========second filter out data count==========" << std::endl;
        purityLog << "#peakValley count: " << filterCounts.peakValley << std::endl;

        purityLog << "#==========wisker filter out data count==========" << std::endl;
        purityLog << "#iteration times: " << iteration_times << std::endl;
        purityLog << "#remove outliers: " << filterCounts.outliers << std::endl;
        purityLog << "#==========purity prediction===========" << std::endl;
        purityLog << "Tumor purity: " << purity << std::endl;
        purityLog << "Data size: " << plotValue.data_size << std::endl;
        purityLog << "Median: " << plotValue.median << std::endl;
        purityLog << "Q1: " << plotValue.q1 << std::endl;
        purityLog << "Q3: " << plotValue.q3 << std::endl;
        purityLog << "IQR: " << plotValue.iqr << std::endl;
        purityLog << "Whiskers: " << plotValue.lowerWhisker << " to " << plotValue.upperWhisker << std::endl;
        purityLog << "Outliers: " << plotValue.outliers << std::endl;

        purityLog.close();
    }catch(const std::exception& e){
        std::cerr << "[ERROR] :" << e.what() << std::endl;
        std::cerr << "[ERROR] : Failed to write purity log" << std::endl;
    }
}

Histogram::Histogram(){
    histogram = std::vector<HistogramData>(1000, HistogramData(0, 0.0));
    total_snp_count = 0;
    max_height = 0;
    data_range = std::make_pair(0, 0);
}

Histogram::~Histogram(){

}

void Histogram::buildHistogram(const std::vector<PurityData>& purityFeatureValueVec){
    try{
        if(purityFeatureValueVec.empty()){
            throw std::invalid_argument("Purity feature value vector is empty");
        }
        total_snp_count = purityFeatureValueVec.size();

        // build histogram
        for(const auto& data : purityFeatureValueVec){
            size_t readCount = data.germlineReadHpCountInNorBam;
            if(readCount >= histogram.size()){
                size_t new_size = histogram.size() * 2;
                if(new_size >= MAX_HISTOGRAM_SIZE){
                    throw std::overflow_error("Read count " + std::to_string(readCount) + " exceeds maximum histogram size");
                }
                histogram.resize(new_size, HistogramData(0, 0.0));
            }
            histogram[readCount].count++;
        }

        if(histogram.empty()){
            throw std::runtime_error("Histogram is empty after building");
        }

    } catch (const std::exception& e) {
        histogram.clear();
        throw std::runtime_error("Failed to build histogram: " + std::string(e.what()));
    }
}

void Histogram::calculateStatistics(){
    // calculate the percentage of data
    try{
        double total_percentage = 0.0;
        bool load_first_non_zero_data = false;

        for(size_t i = 0; i < histogram.size(); i++){
            double tmp = (double)histogram[i].count / (double)total_snp_count;
            total_percentage += tmp;
            
            histogram[i].percentage = total_percentage;

            //update the max height
            if(histogram[i].count > max_height){
                max_height = histogram[i].count;
            }

            //find the first non-zero data
            if(load_first_non_zero_data == false && histogram[i].count > 0){
                data_range.first = i;
                load_first_non_zero_data = true;
            }
            // update the position of the last non-zero value
            if(histogram[i].count > 0){
                data_range.second = i;  
            }
        }

        if(max_height == 0){
            throw std::runtime_error("max_height is 0 in histogram");
        }

        //resize the histogram
        histogram.resize(data_range.second + 1);
        histogram.shrink_to_fit();
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to calculate statistics: " + std::string(e.what()));
    }
}

void Histogram::applyGaussianFilter(double sigma) {
    try{
        if (histogram.empty()) {
            throw std::runtime_error("histogram is empty");
        }

        std::vector<double> kernel = createGaussianKernel(sigma);
        
        if (kernel.empty() || kernel.size() % 2 == 0) {
            throw std::invalid_argument("Invalid kernel: must be non-empty with odd size");
        }

        std::vector<HistogramData> temp = histogram; 
        const size_t half_size = kernel.size() / 2;

        // Apply convolution with edge padding
        for (size_t i = 0; i < histogram.size(); ++i) {
            double smoothed_count = 0.0;

            for (size_t j = 0; j < kernel.size(); ++j) {
                size_t idx = 0;
                if (i + j >= half_size) {
                    idx = i + j - half_size;
                    if (idx >= histogram.size()) {
                        idx = histogram.size() - 1;
                    }
                }

                smoothed_count += temp[idx].count * kernel[j];
            }

            // Check for overflow
            if (!std::isfinite(smoothed_count)) {
                throw std::runtime_error("Invalid smoothed value detected: possible numerical overflow");
            }

            histogram[i].count = smoothed_count;
        }

        // Update the statistics
        calculateStatistics();

    }catch(const std::exception& e){
        throw std::runtime_error("Gaussian filter failed: " + std::string(e.what()));
    }
} 

std::vector<double> Histogram::createGaussianKernel(double sigma) {
    const double MIN_SIGMA = 0.1;
    const double MAX_SIGMA = 2.0;

    if (sigma < MIN_SIGMA || sigma > MAX_SIGMA) {
        throw std::invalid_argument("Sigma must be between 0.1 and 2.0");
    }

    // Calculate kernel size (6 * sigma + 1)
    int kernel_size = static_cast<int>(KERNEL_SIZE_MULTIPLIER * sigma + 1);
    if (kernel_size % 2 == 0) {
        kernel_size += 1;  // Ensure odd kernel size
    }

    std::vector<double> kernel(kernel_size);
    int half_size = kernel_size / 2;

    // Create Gaussian kernel
    double sum = 0.0;
    for (int i = 0; i < kernel_size; ++i) {
        double x = i - half_size;
        kernel[i] = std::exp(-0.5 * (x / sigma) * (x / sigma));
        sum += kernel[i];
    }

    // Normalize kernel
    for (double& k : kernel) {
        k /= sum;
    }

    return kernel;
}

Histogram Histogram::getSmoothedHistogram(double sigma) {
    // create a copy of the current object
    Histogram temp = *this; 
    try{
        // apply gaussian filter
        temp.applyGaussianFilter(sigma);
        return temp;
    }catch(const std::exception& e){
        std::cerr << "[ERROR] Failed to get smoothed histogram: " << e.what() << std::endl;
        std::cerr << "[INFO] Falling back to original histogram" << std::endl;
        return *this;
    }
} 

PeakProcessor::PeakProcessor(){
    mainPeakCount = 0;
    mainPeak = MainPeakInfo();
    saddlePoint = SaddlePointInfo();
}

PeakProcessor::~PeakProcessor(){

}

void PeakProcessor::findPeaks(const std::vector<HistogramData>& histogram, const double min_peak_height){
    try{
        if(histogram.empty()){
            throw std::invalid_argument("Histogram is empty");
        }
        // find peak value and its index
        for(size_t i = 0; i < histogram.size(); i++){

            bool is_peak = false;
            //skip the data which is less than min_peak_count
            if(histogram[i].count < min_peak_height){
                continue;
            }
            //first element
            else if(i == 0 && i != histogram.size() - 1){
                if(histogram[i].count > histogram[i+1].count){
                    is_peak = true;
                }
            }
            //last element
            else if(i == histogram.size() - 1 && i != 0){
                if(histogram[i].count > histogram[i-1].count){
                    is_peak = true;
                }
            }
            else if( histogram[i].count > histogram[i-1].count && histogram[i].count > histogram[i+1].count){
                is_peak = true;
            }

            if(is_peak){
                peaksVec.emplace_back(Peak{i, histogram[i].count});
            }
        }

        //write the peak info to the exec_log
        for(size_t i = 0; i < peaksVec.size(); i++){
            exec_log.push_back("[INFO] Peak " + std::to_string(i) + ": " + std::to_string(peaksVec[i].histo_index) + ", " + std::to_string(peaksVec[i].height));
        }
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to find peaks: " + std::string(e.what()));
    }
}

void PeakProcessor::removeClosePeaks(const size_t minDistance) {
    try{
        if (peaksVec.empty()) {
            throw std::runtime_error("No peaks found in peaksVec");
        }else if(peaksVec.size() >= 2){
            //remove close peak with min_distance
            for(size_t i = 0; i < peaksVec.size() - 1;) {
                if(peaksVec[i+1].histo_index - peaksVec[i].histo_index < minDistance) {
                    if(peaksVec[i].height >= peaksVec[i+1].height) {
                        exec_log.push_back("[INFO] remove the peak " + std::to_string(peaksVec[i+1].histo_index) + "(" + std::to_string(peaksVec[i+1].height) + ")" + " because it is too close to the peak " + std::to_string(peaksVec[i].histo_index) + "(" + std::to_string(peaksVec[i].height) + ")");
                        peaksVec.erase(peaksVec.begin() + i + 1);
                    } else {
                        exec_log.push_back("[INFO] remove the peak " + std::to_string(peaksVec[i].histo_index) + "(" + std::to_string(peaksVec[i].height) + ")" + " because it is too close to the peak " + std::to_string(peaksVec[i+1].histo_index) + "(" + std::to_string(peaksVec[i+1].height) + ")");
                        peaksVec.erase(peaksVec.begin() + i);
                    }
                } else {
                    ++i;
                }
            }
        }
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to remove close peaks: " + std::string(e.what()));
    }
}

void PeakProcessor::determineTrends() {
    try{
        if(peaksVec.empty()){
            throw std::runtime_error("No peaks found in peaksVec");
        }
        //determine the trend of the peak
        if(peaksVec.size() >= 2){
            for(size_t i = 0; i < peaksVec.size() - 1; i++) {
                if(peaksVec[i].height < peaksVec[i+1].height) {
                    peaksVec[i].right_trend = PeakTrend::UP;
                    peaksVec[i+1].left_trend = PeakTrend::UP;
                } else if(peaksVec[i].height > peaksVec[i+1].height) {
                    peaksVec[i].right_trend = PeakTrend::DOWN;
                    peaksVec[i+1].left_trend = PeakTrend::DOWN;
                } else {
                    peaksVec[i].right_trend = PeakTrend::FLAG;
                    peaksVec[i+1].left_trend = PeakTrend::FLAG;
                }
            }
        }
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to determine trends: " + std::string(e.what()));
    }
}

void PeakProcessor::findMainPeakCandidates() {
    try{
        //find the main peak
        if(peaksVec.empty()){
            throw std::runtime_error("No peaks found in peaksVec");
        }else if(peaksVec.size() == 1){
            exec_log.push_back("[INFO] Only one peak found");
            peaksVec[0].is_main_peak = true;
            mainPeakCount = 1;
        }else{
            for(size_t index = 0; index < peaksVec.size(); index++) {
                bool is_main_peak = false;
                //first peak
                if(index == 0) {
                if(peaksVec[index].right_trend == PeakTrend::DOWN) {
                    is_main_peak = true;
                }
                //last peak
                } else if(index == peaksVec.size() - 1) {
                    if(peaksVec[index].left_trend == PeakTrend::UP) {
                        is_main_peak = true;
                    }
                //middle peak
                } else if(peaksVec[index].left_trend == PeakTrend::UP && peaksVec[index].right_trend == PeakTrend::DOWN) {
                    is_main_peak = true;
                }
                if(is_main_peak) {
                    peaksVec[index].is_main_peak = true;
                        mainPeakCount++;
                }
            }
        }
    }catch(const std::exception& e) {
        throw std::runtime_error("Failed to find main peak candidates: " + std::string(e.what()));
    }
}

bool PeakProcessor::findFirstPriorityMainPeak() {

    bool found_first_main_peak = false;
    try{
        std::vector<Peak> mainPeakVec;

        for(auto& peak : peaksVec) {
            if(peak.is_main_peak) {
                mainPeakVec.push_back(peak);
            }
        }

        if(mainPeakVec.empty()) {
            throw std::runtime_error("No main peaks found in peaksVec");
        }

        // only one main peak
        if(mainPeakVec.size() == 1) {
            mainPeak.index =  mainPeakVec[0].histo_index;
            found_first_main_peak = true;
        }else{
            //sort the main peak candidates by height in descending order
            std::sort(mainPeakVec.begin(), mainPeakVec.end(), Peak::compareByHight);
            //get the first and second higher peak
            //seclect higher index peak as the first main peak from the two higher peaks
            if(mainPeakVec[0].histo_index > mainPeakVec[1].histo_index){
                mainPeak.index =  mainPeakVec[0].histo_index;
            }else{
                mainPeak.index =  mainPeakVec[1].histo_index;
            }
            found_first_main_peak = true;
        }
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to find first priority main peak: " + std::string(e.what()));
    }

    return found_first_main_peak;
}

bool PeakProcessor::findSaddlePoint() {
    //find the saddle point
    bool found_saddle_point = false;
    try{
        auto peakVecIter = peaksVec.begin();
        //if the first main peak is the first peak
        if (peakVecIter->histo_index == mainPeak.index){
            saddlePoint.index = -1;
            exec_log.push_back("[INFO] the first main peak is the first peak");
        }else{
            //align the peakVecIter to the first main peak
            while(peakVecIter->histo_index != mainPeak.index){
                if(peakVecIter == peaksVec.end()){
                    throw std::runtime_error("Main peak not found in peaksVec");
                }
                peakVecIter++;
            }

            //align the peakVecIter to the left peak of the first main peak
            peakVecIter--;

            //if the first main peak is the second peak
            if(peakVecIter == peaksVec.begin()){
                saddlePoint.index = peakVecIter->histo_index;
                found_saddle_point = true;
            }else{
                while(peakVecIter != peaksVec.begin()){
                    //record the valley of the peak or the first peak as the saddle point
                    if((peakVecIter->left_trend == PeakTrend::DOWN && peakVecIter->right_trend == PeakTrend::UP)){
                        saddlePoint.index = peakVecIter->histo_index;
                        found_saddle_point = true;
                        break;
                    }
                    peakVecIter--;
                }
                //if no saddle point found, select the first peak as the saddle point
                if(!found_saddle_point){
                    exec_log.push_back("[INFO] no saddle point found, select the first peak as the saddle point");
                    saddlePoint.index = peaksVec.begin()->histo_index;
                    found_saddle_point = true;
                }
            }
        }

    }catch(const std::exception& e){
        throw std::runtime_error("Failed to find saddle point: " + std::string(e.what()));
    }

    return found_saddle_point;
}

bool PeakProcessor::findLowestValley(const std::vector<HistogramData>& histogram, size_t start_index, size_t end_index, Valley& result) {
    // check the index
    if (start_index >= end_index || end_index > histogram.size()) {
        exec_log.push_back("[ERROR] (findLowestValley) index out of range: start: " + std::to_string(start_index) + " end: " + std::to_string(end_index) + " histogram.size(): " + std::to_string(histogram.size()));
        return false;
    }
    
    bool found = false;
    result.height = INT_MAX;
    
    // find the lowest valley
    for(size_t i = start_index + 1; i < end_index - 1; i++) {
        if(histogram[i].count < histogram[i-1].count && 
           histogram[i].count < histogram[i+1].count) {
            if (!found || histogram[i].count < result.height) {
                result.index = i;
                result.height = histogram[i].count;
                result.percentage = histogram[i].percentage;
                found = true;
            }
        }
    }
    
    return found;
}

void PeakProcessor::SetThresholdByValley(const std::vector<HistogramData>& histogram, const double &max_height){
    bool found_first_main_peak = false;
    mainPeak.peak = Peak();
    
    bool found_saddle_point = false;
    saddlePoint.peak = Peak();

    saddlePoint.next_peak = Peak();
    saddlePoint.pre_peak = Peak();

    lowestValley = Valley();
    thresholdPercentage = 0.0;
    threshold = 0;

    //find the first main peak
    found_first_main_peak = findFirstPriorityMainPeak();

    if(found_first_main_peak){
        mainPeak.peak = getPeak(mainPeak.index, 0);
        exec_log.push_back("[INFO] found the first main peak :" + std::to_string(mainPeak.peak.histo_index));
        //find the saddle point
        found_saddle_point = findSaddlePoint();
        if(found_saddle_point){
            saddlePoint.peak = getPeak(saddlePoint.index, 0);
            exec_log.push_back("[INFO] found the saddle point :" + std::to_string(saddlePoint.peak.histo_index));
            exec_log.push_back("[INFO] check the next peak of the saddle point");

            saddlePoint.next_peak = getPeak(saddlePoint.index, 1);
            // find the lowest height valley between the saddle point and its next peak
            bool found_valley = findLowestValley(histogram, saddlePoint.peak.histo_index, saddlePoint.next_peak.histo_index, lowestValley);
            if(found_valley){
                exec_log.push_back("[INFO] find the lowest height valley: " + std::to_string(lowestValley.index) + "(" + std::to_string(lowestValley.percentage) + ")");
                thresholdPercentage = lowestValley.percentage;
                threshold = lowestValley.index;
            }else{
                exec_log.push_back("[INFO] no valley found");
            }


            //threshold >= 0.3, remove the lowest height valley
            if(thresholdPercentage >= THRESHOLD_PERCENTAGE_LIMIT || !found_valley){
                // reset the lowest valley
                lowestValley = Valley();
                thresholdPercentage = 0.0;
                threshold = 0;
                bool found_valley = false;

                exec_log.push_back("[INFO] threshold >= " + std::to_string(THRESHOLD_PERCENTAGE_LIMIT) + "%, reset threshold to " + std::to_string(threshold) + "(" + std::to_string(thresholdPercentage) + ")");
                exec_log.push_back("[INFO] check the pre peak of the saddle point");
                
                // saddle point is not the first peak
                if(saddlePoint.peak.histo_index != peaksVec[0].histo_index){
                    saddlePoint.pre_peak = getPeak(saddlePoint.peak.histo_index, -1);
                    exec_log.push_back("[INFO] saddle point have a pre peak " + std::to_string(saddlePoint.pre_peak.histo_index) + "->" + std::to_string(saddlePoint.peak.histo_index));
                    //find the lowest height valley between the saddle point and its pre peak
                    found_valley = findLowestValley(histogram, saddlePoint.pre_peak.histo_index, saddlePoint.peak.histo_index, lowestValley);
                    if(found_valley){
                        exec_log.push_back("[INFO] find the lowest height valley : " + std::to_string(lowestValley.index) + "(" + std::to_string(lowestValley.percentage) + ")");
                        thresholdPercentage = lowestValley.percentage;
                        threshold = lowestValley.index;
                    }else{
                        exec_log.push_back("[INFO] no valley found between the saddle point and its pre peak");
                    }
                }else{
                    exec_log.push_back("[INFO] no pre peak found");
                }
            }
        }else{
            exec_log.push_back("[INFO] no saddle point found");
        }
    }
    //check valley height 
    //valley may not be able to separate the two distributions
    if(lowestValley.height > max_height * 0.7){
        exec_log.push_back("[INFO] valley height is too high, set the threshold to 0: valley height: " + std::to_string(lowestValley.height) + " max height: " + std::to_string(max_height));
        thresholdPercentage = 0.0;
        threshold = 0;
    }

    // Final check threshold
    if(thresholdPercentage >= THRESHOLD_PERCENTAGE_LIMIT){
        exec_log.push_back("[INFO] Final threshold over " + std::to_string(THRESHOLD_PERCENTAGE_LIMIT) + "%, set to 0: "+ std::to_string(threshold) + "(" + std::to_string(thresholdPercentage) + ")");
        thresholdPercentage = 0.0;
        threshold = 0;
    }
    exec_log.push_back("[INFO] Final threshold: " + std::to_string(threshold) + "(" + std::to_string(thresholdPercentage) + ")");

}

Peak PeakProcessor::getPeak(size_t histo_index, int offset){
    try{
        for(size_t i = 0; i < peaksVec.size(); i++){
            if(peaksVec[i].histo_index == histo_index){
                if(i + offset < peaksVec.size() || i + offset >= 0){
                    return peaksVec[i + offset];
                }else{
                    throw std::runtime_error("Peak index out of range: peaksVec.size(): " + std::to_string(peaksVec.size()) + " index: " + std::to_string(i) + " offset: " + std::to_string(offset));
                }
            }
        }
        throw std::runtime_error("Peak not found: histo_index: " + std::to_string(histo_index) + " offset: " + std::to_string(offset));
    }catch(const std::exception& e){
        throw std::runtime_error("Failed to get peak: " + std::string(e.what()));
    }
}

int PeakProcessor::getThreshold(){
    return threshold;
}

std::string PeakProcessor::transformTrend(const PeakTrend &trend) {
    switch(trend) {
        case PeakTrend::NONE: return "NONE";
        case PeakTrend::UP: return "UP";
        case PeakTrend::DOWN: return "DOWN";
        case PeakTrend::FLAG: return "FLAG";
        default: return "UNKNOWN";
    }
}

void PeakProcessor::writePeakValleyLog(
    const HaplotagParameters &params,
    const std::vector<HistogramData>& histogram,
    const std::vector<HistogramData>& smoothed_histogram,
    size_t &total_snp_count,
    const std::pair<size_t, size_t>& data_range,
    double &max_height,
    const double &MIN_PEAK_RATIO,
    double &peak_threshold,
    double &sigma) {

    std::string postfix = "_germlineReadHpCountInNorBam_histogram.out";
    std::ofstream histogramFile(params.resultPrefix + postfix);
    if(!histogramFile.is_open()){
        std::cerr << "[WARNING] Failed to open histogram log file: " << params.resultPrefix + postfix << std::endl;        
        return;
    }

    histogramFile << "#total snp count: " << total_snp_count << std::endl;
    histogramFile << "#data range: " << data_range.first << " to " << data_range.second << std::endl;
    histogramFile << "#max height: " << max_height << std::endl;
    histogramFile << "#min peak ratio: " << MIN_PEAK_RATIO << std::endl;
    histogramFile << "#peak threshold: " << peak_threshold << std::endl;
    histogramFile << "#gaussian filter sigma: " << sigma << std::endl;

    //write exec log to file
    histogramFile << "#========Execution Log==========" << std::endl;
    for(const auto& log : exec_log){
        histogramFile << "#" << log << std::endl;
    }
    histogramFile << "\nindex (germline Hp read count in normal bam), height (snp count), percentage\n" << std::endl;
    //write smoothed histogram 
    histogramFile << "#Smoothed Histogram Start" << std::endl;
    for(size_t i = 0; i < smoothed_histogram.size(); i++){
        histogramFile << std::fixed << std::setprecision(2) << i << "\t" << smoothed_histogram[i].count << "\t" << smoothed_histogram[i].percentage << std::endl;
    }
    histogramFile << "#Smoothed Histogram End\n" << std::endl;

    //write histogram 
    histogramFile << "#Histogram Start" << std::endl;
    for(size_t i = 0; i < histogram.size(); i++){
        histogramFile << i << "\t" << histogram[i].count << "\t" << histogram[i].percentage << std::endl;
    }
    histogramFile << "#Histogram End\n" << std::endl;

    histogramFile << "\n#==========Peak Trend Analysis==========" << std::endl;
    histogramFile << "#peak count: " << peaksVec.size() << std::endl;
    histogramFile << "#Peak " << "\t" 
                  << std::left << std::setw(10) << "Position" << "\t"
                  << std::left << std::setw(10) << "Height" << "\t"
                  << std::left << std::setw(10) << "Left_Trend" << "\t"
                  << std::left << std::setw(10) << "Right_Trend" << std::endl;    
    for(size_t i = 0; i < peaksVec.size(); i++){
        std::string left_trend = transformTrend(peaksVec[i].left_trend);
        std::string right_trend = transformTrend(peaksVec[i].right_trend);
        histogramFile << i+1 << "\t" 
                    << std::left << std::setw(10) << peaksVec[i].histo_index << "\t" 
                    << std::left << std::setw(10) << peaksVec[i].height << "\t" 
                    << std::left << std::setw(10) << left_trend << "\t" 
                    << std::left << std::setw(10) << right_trend << std::endl;
    }
    
    histogramFile << "\n#==========Main Peak Analysis==========" << std::endl;
    if(mainPeakCount == 0){
        histogramFile << "#main peak count: 0" << std::endl;
    }else{
        histogramFile << "#main peak count: " << mainPeakCount << std::endl;
        histogramFile << "#Peak " << "\t" 
                      << std::left << std::setw(10) << "Position" << "\t"
                      << std::left << std::setw(10) << "Height" << "\t"
                      << std::left << std::setw(10) << "Left_Trend" << "\t"
                      << std::left << std::setw(10) << "Right_Trend" << std::endl;
        auto peakVecIter = peaksVec.begin();
        int i = 1;
        while(peakVecIter != peaksVec.end()){
            if(peakVecIter->is_main_peak){
                histogramFile << i << "\t" 
                              << std::left << std::setw(10) << peakVecIter->histo_index << "\t" 
                              << std::left << std::setw(10) << peakVecIter->height << "\t" 
                              << std::left << std::setw(10) << transformTrend(peakVecIter->left_trend) << "\t" 
                              << std::left << std::setw(10) << transformTrend(peakVecIter->right_trend) << std::endl;
            }
            i++;
            peakVecIter++;
        } 
    }

    histogramFile << "\n#==========Selected Peak & Valley==========" << std::endl;
    histogramFile << "#first main peak       : " << mainPeak.peak.histo_index << "\t" << mainPeak.peak.height << "\t" << transformTrend(mainPeak.peak.left_trend) << "\t" << transformTrend(mainPeak.peak.right_trend) << std::endl;
    histogramFile << "#saddle point          : " << saddlePoint.peak.histo_index << "\t" << saddlePoint.peak.height << "\t" << transformTrend(saddlePoint.peak.left_trend) << "\t" << transformTrend(saddlePoint.peak.right_trend) << std::endl;
    histogramFile << "#saddle point next peak: " << saddlePoint.next_peak.histo_index << "\t" << saddlePoint.next_peak.height << "\t" << transformTrend(saddlePoint.next_peak.left_trend) << "\t" << transformTrend(saddlePoint.next_peak.right_trend) << std::endl;
    histogramFile << "#saddle point pre peak : " << saddlePoint.pre_peak.histo_index << "\t" << saddlePoint.pre_peak.height << "\t" << transformTrend(saddlePoint.pre_peak.left_trend) << "\t" << transformTrend(saddlePoint.pre_peak.right_trend) << std::endl;
    histogramFile << "#lowest height valley  : " << lowestValley.index << "\t" << lowestValley.percentage << std::endl;
    histogramFile << "#threshold percentage: " << thresholdPercentage << std::endl;
    histogramFile << "#threshold: " << threshold << std::endl;
    histogramFile.close();  
}