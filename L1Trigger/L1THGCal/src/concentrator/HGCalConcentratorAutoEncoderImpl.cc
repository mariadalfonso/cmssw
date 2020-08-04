#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorAutoEncoderImpl.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include <iomanip>

// Following example of implementing graphloading from here:
// https://gitlab.cern.ch/mrieger/CMSSW-TensorFlowExamples/-/blob/master/GraphLoading/

HGCalConcentratorAutoEncoderImpl::HGCalConcentratorAutoEncoderImpl(const edm::ParameterSet& conf)
    : cellRemap_(conf.getParameter<std::vector<int>>("cellRemap")),
      bitsPerInput_(conf.getParameter<int>("nBitsPerInput")),
      maxBitsPerOutput_(conf.getParameter<int>("maxBitsPerOutput")),
      outputBitsPerLink_(conf.getParameter<std::vector<int>>("bitsPerLink")),
      modelFilePaths_(conf.getParameter<std::vector<edm::ParameterSet>>("modelFiles")),
      linkToGraphMap_(conf.getParameter<std::vector<unsigned int>>("linkToGraphMap")),
      zeroSuppresionThreshold_(conf.getParameter<double>("zeroSuppresionThreshold")),
      saveEncodedValues_(conf.getParameter<bool>("saveEncodedValues")),
      preserveModuleSum_(conf.getParameter<bool>("preserveModuleSum")) {
  //construct inverse array, to get U/V for a particular ae output position
  for (unsigned i = 0; i < cellRemap_.size(); i++) {
    if (cellRemap_[i] > -1) {
      if (cellRemap_[i] > nTriggerCells_ - 1) {
        throw cms::Exception("BadInitialization")
            << "cellRemap value " << cellRemap_[i] << " is larger than the number of trigger cells " << nTriggerCells_;
      }
      ae_outputCellU_[cellRemap_[i]] = int(i / cellRemapRowOffset_);
      ae_outputCellV_[cellRemap_[i]] = i % cellRemapRowOffset_;
    }
  }

  tensorflow::setLogging("0");

  for (const auto& modelFilePset : modelFilePaths_) {
    std::string encoderPath = modelFilePset.getParameter<edm::FileInPath>("encoderModelFile").fullPath();
    std::string decoderPath = modelFilePset.getParameter<edm::FileInPath>("decoderModelFile").fullPath();

    graphDef_encoder_ = std::unique_ptr<tensorflow::GraphDef>{tensorflow::loadGraphDef(encoderPath)};

    // create a new session and add the graphDef
    session_encoder_.push_back(
        std::unique_ptr<tensorflow::Session>{tensorflow::createSession(graphDef_encoder_.get())});

    graphDef_decoder_ = std::unique_ptr<tensorflow::GraphDef>{tensorflow::loadGraphDef(decoderPath)};

    // create a new session and add the graphDef
    session_decoder_.push_back(
        std::unique_ptr<tensorflow::Session>{tensorflow::createSession(graphDef_decoder_.get())});

    //extract encoder tenser names from first graph, check that rest of the names are consistent
    if (modelFilePset == modelFilePaths_.front()) {
      inputTensorName_encoder_ = graphDef_encoder_.get()->node(0).name();
      outputTensorName_encoder_ = graphDef_encoder_.get()->node(graphDef_encoder_.get()->node_size() - 1).name();
      inputTensorName_decoder_ = graphDef_decoder_.get()->node(0).name();
      outputTensorName_decoder_ = graphDef_decoder_.get()->node(graphDef_decoder_.get()->node_size() - 1).name();
    } else {
      if (inputTensorName_encoder_ != graphDef_encoder_.get()->node(0).name()) {
        throw cms::Exception("BadInitialization") << "provided list of encoder graphs have different input nodes";
      }
      if (outputTensorName_encoder_ != graphDef_encoder_.get()->node(graphDef_encoder_.get()->node_size() - 1).name()) {
        throw cms::Exception("BadInitialization") << "provided list of encoder graphs have different output nodes";
      }
      if (inputTensorName_decoder_ != graphDef_decoder_.get()->node(0).name()) {
        throw cms::Exception("BadInitialization") << "provided list of decoder graphs have different input nodes";
      }
      if (outputTensorName_decoder_ != graphDef_decoder_.get()->node(graphDef_decoder_.get()->node_size() - 1).name()) {
        throw cms::Exception("BadInitialization") << "provided list of decoder graphs have different output nodes";
      }
    }
  }

  // check that the appropriate number of links have been specified
  if (linkToGraphMap_.size() <= maxNumberOfLinks_) {
    throw cms::Exception("BadInitialization")
        << "Autoencoder graph number must be specified for all link allocation possibilities. Only "
        << linkToGraphMap_.size() << " values specified while " << maxNumberOfLinks_ << "links are possible";
  }

  // check that all graph indices specified exist in the model file lists
  for (const auto& graphNumber : linkToGraphMap_) {
    if (graphNumber >= modelFilePaths_.size()) {
      throw cms::Exception("BadInitialization")
          << "Autoencoder graph number  " << graphNumber << " is larger than the size of the provided list of graphs "
          << modelFilePaths_.size();
    }
  }
}

void HGCalConcentratorAutoEncoderImpl::select(unsigned nLinks,
                                              const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
                                              std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput,
                                              std::vector<l1t::HGCalConcentratorData>& ae_encodedLayer_Output) {
  std::array<double, nTriggerCells_> mipPt;
  std::array<double, nTriggerCells_> uncompressedCharge;
  std::array<double, nTriggerCells_> compressedCharge;
  std::array<double, nTriggerCells_> ae_inputArray;
  std::array<double, nTriggerCells_> ae_outputArray;

  //reset inputs to 0 to account for zero suppressed trigger cells
  mipPt.fill(0);
  uncompressedCharge.fill(0);
  compressedCharge.fill(0);

  double modSum = 0;

  int bitsPerOutput = outputBitsPerLink_.at(nLinks);

  // largest expected input and output values, used for bit truncation
  // values of -1 for the number of bits used to keep full precision, in which case the MaxIntSize variables are not used
  double inputMaxIntSize = 1;
  if (bitsPerInput_ > 0)
    inputMaxIntSize = 1 << bitsPerInput_;
  double outputMaxIntSize = 1;
  if (bitsPerOutput > 0)
    outputMaxIntSize = 1 << bitsPerOutput;
  double outputMaxIntSizeGlobal = 1;
  if (maxBitsPerOutput_ > 0)
    outputMaxIntSizeGlobal = 1 << maxBitsPerOutput_;

  for (const auto& trigCell : trigCellVecInput) {
    if (triggerTools_.isScintillator(trigCell.detId()))
      return;  //currently, only silicon modules are setup to work (mapping of scinillators would be different, and needs to be investigated)

    HGCalTriggerDetId id(trigCell.detId());
    uint cellu = id.triggerCellU();
    uint cellv = id.triggerCellV();
    int inputIndex = cellRemap_.at(cellu * cellRemapRowOffset_ + cellv);
    if (inputIndex < 0) {
      throw cms::Exception("BadInitialization")
          << "Invalid index provided for trigger cell u=" << cellu << " v=" << cellv << " in cellRemap["
          << (cellu * cellRemapRowOffset_ + cellv) << "]";
    }

    mipPt[inputIndex] = trigCell.mipPt();
    uncompressedCharge[inputIndex] = trigCell.uncompressedCharge();
    compressedCharge[inputIndex] = trigCell.compressedCharge();

    modSum += trigCell.mipPt();
  }

  if (modSum > 0) {
    //normalize inputs to module sum
    for (int i = 0; i < nTriggerCells_; i++) {
      ae_inputArray[i] = mipPt[i] / modSum;
      //round to precision of input, if bitsPerInput_ is -1 keep full precision
      if (bitsPerInput_ > 0) {
        ae_inputArray[i] = std::round(ae_inputArray[i] * inputMaxIntSize) / inputMaxIntSize;
      }
    }
  }

  tensorflow::Tensor encoder_input(tensorflow::DT_FLOAT,
                                   {encoderShape_0_, encoderShape_1_, encoderShape_2_, encoderShape_3_});
  float* d = encoder_input.flat<float>().data();
  for (int i = 0; i < nTriggerCells_; i++, d++) {
    *d = ae_inputArray[i];
  }

  int graphIndex = linkToGraphMap_.at(nLinks);

  std::vector<tensorflow::Tensor> encoder_outputs;
  tensorflow::run(session_encoder_.at(graphIndex).get(),
                  {{inputTensorName_encoder_, encoder_input}},
                  {outputTensorName_encoder_},
                  &encoder_outputs);

  if (encoder_outputs.empty()) {
    throw cms::Exception("BadInitialization") << "Autoencoder graph returning empty output vector";
  }

  d = encoder_outputs[0].flat<float>().data();
  for (int i = 0; i < encoder_outputs[0].NumElements(); i++, d++) {
    ae_encodedLayer_[i] = *d;
    //truncate the encoded layer bits
    if (bitsPerOutput > 0) {
      ae_encodedLayer_[i] = std::round(ae_encodedLayer_[i] * outputMaxIntSize) / outputMaxIntSize;
    }
  }

  tensorflow::Tensor decoder_input(tensorflow::DT_FLOAT, {decoderShape_0_, decoderShape_1_});
  d = decoder_input.flat<float>().data();
  for (int i = 0; i < nEncodedLayerNodes_; i++, d++) {
    *d = ae_encodedLayer_[i];
  }

  std::vector<tensorflow::Tensor> decoder_outputs;
  tensorflow::run(session_decoder_.at(graphIndex).get(),
                  {{inputTensorName_decoder_, decoder_input}},
                  {outputTensorName_decoder_},
                  &decoder_outputs);

  double outputSum = 1.;
  d = decoder_outputs[0].flat<float>().data();
  if (preserveModuleSum_) {
    outputSum = 0.;
    for (int i = 0; i < decoder_outputs[0].NumElements(); i++, d++) {
      outputSum += *d;
    }
  }

  d = decoder_outputs[0].flat<float>().data();
  for (int i = 0; i < decoder_outputs[0].NumElements(); i++, d++) {
    ae_outputArray[i] = *d / outputSum;
  }

  // Add data back into trigger cells
  if (modSum > 0) {
    //get detID for everything but cell, take first entry detID and subtract off cellU and cellV contribution
    HGCalTriggerDetId id(trigCellVecInput.at(0).detId());
    int subdet = id.subdet();
    int zp = id.zside();
    int type = id.type();
    int layer = id.layer();
    int waferU = id.waferU();
    int waferV = id.waferV();
    int cellU = id.triggerCellU();
    int cellV = id.triggerCellV();

    //use first TC to find mipPt conversions to Et and ADC
    float mipPtToEt_conv = trigCellVecInput[0].et() / trigCellVecInput[0].mipPt();
    float mipToADC_conv = trigCellVecInput[0].hwPt() / (trigCellVecInput[0].mipPt() * cosh(trigCellVecInput[0].eta()));

    for (int i = 0; i < nTriggerCells_; i++) {
      if (ae_outputArray[i] > 0) {
        cellU = ae_outputCellU_[i];
        cellV = ae_outputCellV_[i];

        HGCalTriggerDetId id(subdet, zp, type, layer, waferU, waferV, cellU, cellV);

        if (!triggerTools_.getTriggerGeometry()->validTriggerCell(id))
          continue;

        GlobalPoint point = triggerTools_.getTCPosition(id);

        double mipPt = ae_outputArray[i] * modSum;
        double adc = mipPt * cosh(point.eta()) * mipToADC_conv;
        double et = mipPt * mipPtToEt_conv;

        if (mipPt < zeroSuppresionThreshold_)
          continue;

        l1t::HGCalTriggerCell triggerCell(reco::LeafCandidate::LorentzVector(), adc, 0, 0, 0, id);
        //Keep the pre-autoencoder charge for this cell
        triggerCell.setUncompressedCharge(uncompressedCharge[i]);
        triggerCell.setCompressedCharge(compressedCharge[i]);
        triggerCell.setMipPt(mipPt);

        math::PtEtaPhiMLorentzVector p4(et, point.eta(), point.phi(), 0.);

        triggerCell.setP4(p4);
        triggerCell.setPosition(point);

        trigCellVecOutput.push_back(triggerCell);
      }
    }

    if (saveEncodedValues_) {
      id = HGCalTriggerDetId(subdet, zp, type, layer, waferU, waferV, 0, 0);
      for (int i = 0; i < nEncodedLayerNodes_; i++) {
        l1t::HGCalConcentratorData encodedLayerData(ae_encodedLayer_[i] * outputMaxIntSizeGlobal, i, id);
        ae_encodedLayer_Output.push_back(encodedLayerData);
      }
    }
  }
}
