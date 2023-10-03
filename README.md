
# Osteoporotic Precise Screening using Chest radiography and Artificial neural network (OPSCAN)

The Osteoporotic Precise Screening using Chest radiography and Artificial neural network (OPSCAN) program aims to identify individuals at high risk for osteoporosis by analyzing chest X-rays performed for non-specific indications using a deep learning model. The goal is to prompt these individuals to undergo DXA testing for confirmation of osteoporosis and early initiation of treatment, thereby preventing subsequent fracture events. The program consists of two parts. The first part involves training the deep learning model using retrospective data. We propose matching chest X-rays with DXA-based T-scores obtained within 180 days for data analysis. Due to privacy concerns, we are unable to provide real patient data and must use your own data for training purposes. However, the syntax for model training is available. The second part involves conducting a randomized controlled trial (registered at clinicaltrial.gov: [NCT05721157](https://classic.clinicaltrials.gov/ct2/show/NCT05721157)). The main research hypothesis is that the OPSCAN program can reduce the problem of underdiagnosis of osteoporosis. Therefore, the primary endpoint of the study is the identification of new-onset osteoporosis patients (confirmed by DXA). This repository is carried out using the software environment R version 3.4.4 with package MXNet version 1.3.0.

The overall file structure is as follows:

```shell
OPSCAN
├── result
│   ├── ...
├── code
│   ├── model traning
│   │   ├── ...
│   ├── RCT analysis
│   │   ├── ...
├── data
│   ├── model traning
│   │   ├── label.csv
│   │   ├── png
│   │   │   ├── ...
│   ├── RCT analysis
│   │   ├── RCT data.RData
```

The relevant syntax for model training and data analysis is provided below. 

## Model training

In the [data/model training](https://github.com/xup6fup/OPSCAN/blob/master/data/model%20training) directory, you will find a file named [label.csv](https://github.com/xup6fup/OPSCAN/blob/master/data/model%20training/label.csv). This dataset is a subset derived from the [ChestX-ray14](https://nihcc.app.box.com/v/ChestXray-NIHCC) database, encompassing chest X-rays in PA view from 100 patients. Within the file, the "NIH_Image_Index" can be used for cross-referencing with the original dataset. Data preprocessing can be conducted using the [D01. pre-processing.R](https://github.com/xup6fup/OPSCAN/blob/master/code/model%20training/D01.%20pre-processing.R), which resizes the CXRs while maintaining their aspect ratio, ensuring that the shorter side is adjusted to 256 pixels. After processing, the results can be found in the "data/model training/RData" directory, which can be used for subsequent deep learning model training. It's important to note that the T-score in this dataset was predicted using the final version of our AI-CXR model. Osteoporosis (designated as the "OP" column) is defined as a T-score <= -2.5. In future applications of these scripts, one can replace the images and labels with their own data for model training.

<img src="data/model training/png/00000001_000.png">

When training models using Mxnet framework, it is essential to first script an iterator. During our model training, we opted to incorporate random cropping of a 224 × 224-pixel region as input, combined with a 50% likelihood of implementing a random horizontal flip. At the inference phase, we utilized a 10-crop evaluation approach, resulting in ten distinct probabilities for each CXR. The ultimate prediction was determined by averaging these ten probabilities. Related code can be found in [M01. cxr_process_core.R](https://github.com/xup6fup/OPSCAN/blob/master/code/model%20training/M01.%20cxr_process_core.R)

Next, it is necessary for us to design the architecture of the deep learning model. You can select a model from the model zoo we provide. In this instance, we utilized the resnet-18 due to its minimal size and rapid execution speed. However, you are free to substitute it with any other model of your preference. You should download a model and put it into the "model" directory. These models have already been pre-trained on ImageNet, thus, with the appropriate parameters, they can be directly employed for image classification tasks encompassing 1,000 categories. Related code can be found in [M02. architecture.R](https://github.com/xup6fup/OPSCAN/blob/master/code/model%20training/M02.%20architecture.R)

Subsequently, we can initiate the training process. By utilizing the command [M03. run.R](https://github.com/xup6fup/OPSCAN/blob/master/code/model%20training/M03.%20run.R), the process can be directly executed. To optimize our deep learning model, we use cross-entropy loss. We incorporated an oversampling technique, which relied on class weights derived from the prevalence of each class within the development set. The Adam optimizer, with its standard parameters, was employed for fine-tuning all network parameters. We operated with a batch size of 32 and initiated with a learning rate of 0.001. This learning rate would be decreased by a factor of 10 whenever the AUC on the tuning set plateaued after an epoch. As a countermeasure against overfitting, early stopping was employed. Specifically, the network was saved after each epoch, and among these, the network model with the least loss on the tuning subset was selected. We commenced our training with a learning rate of 0.001 (referred to as stage 1). If the loss on the tuning set ceased to decrease, we transitioned to a learning rate of 0.0001 (stage 2), followed by a learning rate of 0.00001 (stage 3). The sole regularization approach to counteract overfitting was the incorporation of a weight decay set at 10^-4.

Upon completion of the training, you should be able to locate two files, "OPSCAN-0000.params" and "OPSCAN-symbol.json", within the "model" directory. To analyze new chest X-ray images, simply execute the code [P01. predicting.R](https://github.com/xup6fup/OPSCAN/blob/master/code/model%20training/P01.%20predicting.R). For illustrative purposes, we use the same data for both training and validation in this example. However, it's imperative to note that this is merely a demonstration. One should avoid intermixing training and testing datasets in genuine applications.

## Analysis for the randomized controlled trial

In our randomized controlled trial (RCT), a total of 4,912 individuals were identified by AI-CXR as being at high risk for osteoporosis. These individuals were randomly allocated into either the screening group or the control group. Within the screening group, research assistants relayed the results of the opportunistic AI-enabled CXR to the patients via phone calls. As part of the OPSCAN study, these patients were offered fully reimbursed DXA examinations, encompassing both the lumbar spine and hip, for a half of year. The de-identified data is stored in [RCT data.RData](https://github.com/xup6fup/OPSCAN/blob/master/data/RCT%20analysis/RCT%20data.RData), and the respective analysis scripts are stored in [code/RCT analysis](https://github.com/xup6fup/OPSCAN/blob/master/code/RCT%20analysis) directory.

<img src="result/Fig 02.png">

Above figure provides an overview of the response of patients in the screening group to the OPSCAN program. Within the entire screening group, 12.8% (315 individuals) availed themselves of the free DXA examination facilitated by the OPSCAN program, while an additional 1.8% (43 individuals) underwent DXA examinations arranged by their primary care physicians. The remaining 2,098 individuals, upon being informed of their high-risk status, declined OPSCAN supported examination. Among the 315 patients who underwent the free DXA examination, 75.2% (237 individuals) were diagnosed with osteoporosis, while 23.8% had osteopenia. Only 1.0% of patients had a T score > -1, meaning the positive predictive value of AI-CXR was 75.2%. Of the 237 individuals diagnosed with osteoporosis, 31.6% ultimately opted for anti-osteoporotic treatment following discussions with their physicians.

<img src="result/Fig 03.png">

Above figure compares overall osteoporosis-related outcomes between the OPSCAN program and usual care, analyzed on an intention-to-treat basis. Regarding the primary endpoint, the screening group had 272 individuals (11.1%) newly diagnosed with osteoporosis, while the control group had only 27 individuals (1.1%). This indicates that the OPSCAN program is 11.20 times more effective (95% CI: 7.51-16.71) at identifying previously undiagnosed individuals with osteoporosis than usual care. The higher proportion of high-risk individuals undergoing DXA examinations, facilitated by the free DXA examinations provided through the OPSCAN program, contributed to an 11.47-fold increase (95% CI: 8.10-16.24) in DXA utilization. Additionally, a greater proportion of individuals in the screening group initiated anti-osteoporotic treatment (OR: 5.27; 95% CI: 3.17-8.76) than in the control group.

<img src="result/Fig 04.png">

In above figure panel A, the benefits of the OPSCAN program versus usual care were compared on an intention-to-treat basis in the two sex/age subgroups. Notably, the OPSCAN program demonstrated a significantly greater advantage over usual care in identifying new-onset osteoporosis in the not meeting ISCD indication subgroup (OR: 23.22, 95% CI: 10.15-53.10) than it did in the meeting ISCD indication subgroup (OR: 7.97, 95% CI: 5.03-12.64, p for interaction = 0.027). Similar trends were observed in DXA examination (p for interaction = 0.040) and anti-osteoporotic treatment (p for interaction = 0.078). Panel B illustrates the changes in major osteoporotic risk and hip fracture risk among patients who received DXA examinations by OPSCAN, with and without FRAX score calculation. Prior to the DXA examination, it is evident that the not meeting ISCD indication subgroup had lower risks than the meeting ISCD indication subgroup. However, after undergoing DXA, the not meeting ISCD indication subgroup experienced a substantial increase in risk, with an average increase of 2.10-fold and 7.04-fold in major osteoporotic risk and hip fracture risk, respectively, which was significantly higher than the increase observed in the meeting ISCD indication subgroup (p < 0.001).

<img src="result/Fig 05.png">

Lastly, we assessed the impact of AI-CXR risk stratification on osteoporosis outcomes. It's imperative to emphasize that, given the low-risk group received the same "usual care" as the control group, we excluded the results from the screening group in this context, ensuring a direct comparison only between the low-risk and control groups. This analysis parallels the findings presented in Figure 2.

## Related publications

The paper of OPSCAN is currently under review. The following reference provideS details on similar model training approaches.

