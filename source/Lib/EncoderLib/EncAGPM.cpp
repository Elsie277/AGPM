#include "CommonLib/UnitTools.h"
#include "EncAGPM.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
#if AGPM_ENABLED

void PU::getAGPMMergeCandidates(const PredictionUnit &pu, MergeCtx &mrgCtx)
{
	MergeCtx tmpMergeCtx;

	const Slice &slice = *pu.cs->slice;
	const uint32_t maxNumMergeCand = pu.cs->sps->getMaxNumMergeCand(); //����Merge��ѡ����

	mrgCtx.numValidMergeCand = 0;

	//AGPM��ѡ�б�����Ϊ6�����г�ʼ��
	for (int32_t i = 0; i < ADAPTIVE_PARTITION_MAX_NUM_UNI_CANDS; i++)
	{
		mrgCtx.bcwIdx[i] = BCW_DEFAULT; //˫��Ԥ��Ȩ��
		mrgCtx.interDirNeighbours[i] = 0; //���ڿ��ѡ��Ԥ�ⷽ�򣺵������˫��1��ʾǰ��2��ʾ����3��ʾ˫��
		mrgCtx.mvFieldNeighbours[(i << 1)].refIdx = NOT_VALID; //��ѡ��ǰ��ο�����
		mrgCtx.mvFieldNeighbours[(i << 1) + 1].refIdx = NOT_VALID; //����ο�����
		mrgCtx.mvFieldNeighbours[(i << 1)].mv = Mv(); //��ѡ��ǰ��MV
		mrgCtx.mvFieldNeighbours[(i << 1) + 1].mv = Mv(); //����MV

		mrgCtx.useAltHpelIf[i] = false; 
	}
	//������ͨMerge�б�
	PU::getInterMergeCandidates(pu, tmpMergeCtx, 0);
	
	//����ÿһ����ͨMerge�б��ĺ�ѡ��ѡ�����ʵĺ�ѡ���ӽ�AGPM�ĵ���Ԥ���б�
	for (int32_t i = 0; i < maxNumMergeCand; i++)
	{	
		int parity = i & 1; //��ż�ԣ�1���棬0��ż

		//ѡ����ͨMerge�б���ż����ѡ��ǰ��MV or ������ѡ�ĺ���MV
		if (tmpMergeCtx.interDirNeighbours[i] & (0x01 + parity))
		{
			mrgCtx.interDirNeighbours[mrgCtx.numValidMergeCand] = 1 + parity; //��������merge��ѡ��Ӧ��Ԥ�ⷽ��(ǰ ��)
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + !parity].mv = Mv(0, 0); //ż�ĺ���MVΪ(0,0),���ǰ��Ϊ(0,0)
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + parity].mv = tmpMergeCtx.mvFieldNeighbours[(i << 1) + parity].mv;
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + !parity].refIdx = -1; //�ο�֡����
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + parity].refIdx = tmpMergeCtx.mvFieldNeighbours[(i << 1) + parity].refIdx;
			
			mrgCtx.numValidMergeCand++;

			//i��mrgCtx.numValidMergeCand��ͬ������������AGPM��ѡ�б�������ѭ��
			if (mrgCtx.numValidMergeCand == ADAPTIVE_PARTITION_MAX_NUM_UNI_CANDS)
			{
				return;
			}
			continue;
		}

		//����ͨMerge��ż��ʱû��ǰ����ȡ����MV�����������޺����ѡ����ȡǰ���ѡ
		if (tmpMergeCtx.interDirNeighbours[i] & (0x02 - parity))
		{
			mrgCtx.interDirNeighbours[mrgCtx.numValidMergeCand] = 2 - parity; //��ѡ��Ԥ�ⷽ��
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + !parity].mv = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].mv;
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + parity].mv = Mv(0, 0);
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + !parity].refIdx = tmpMergeCtx.mvFieldNeighbours[(i << 1) + !parity].refIdx;
			mrgCtx.mvFieldNeighbours[(mrgCtx.numValidMergeCand << 1) + parity].refIdx = -1;

			mrgCtx.numValidMergeCand++;

			if (mrgCtx.numValidMergeCand == ADAPTIVE_PARTITION_MAX_NUM_UNI_CANDS)
			{
				return;
			}
		}
	}
}

void EncCu::xCheckRDCostMergeAGPM2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner,
	const EncTestMode &encTestMode)
{
	const Slice &slice = *tempCS->slice;
	const SPS   &sps = *tempCS->sps;
	CHECK(slice.getSliceType() != B_SLICE, "Adaptive geo mode is only applied to B-slices");

	//��ʼ�����ݽṹ
	tempCS->initStructData(encTestMode.qp);

	bool AdpPartitionHasNoResidual[ADAPTIVE_PARTITION_MAX_NUM_CANDS]; //��Ϻ�ѡ�Ƿ��вв�,trianglecandHasNoResidual[TRIANGLE_MAX_NUM_CANDS] һ�ֻ��ַ�ʽ*6*5
	for (int mergeCand = 0; mergeCand < ADAPTIVE_PARTITION_MAX_NUM_CANDS; mergeCand++)
	{
		AdpPartitionHasNoResidual[mergeCand] = false; //��ʼ��Ϊfalse
	}

	bool bestIsSkip = false;

	uint8_t												  numAdpPartitionCand = ADAPTIVE_PARTITION_MAX_NUM_CANDS; // AGPM��ѡ������30
	uint8_t											      adpPartitionNumMrgSATDCand = ADAPTIVE_PARTITION_MAX_NUM_SATD_CANDS; //����SATD�Ƚϵĺ�ѡ�����Ϊ3
	PelUnitBuf											  adpPartitionBuffer[ADAPTIVE_PARTITION_MAX_NUM_UNI_CANDS]; //���浥Ԫ�����Ϊ6����Ԥ��ĵ���Merge�б�����Ϊ6
	PelUnitBuf											  adpPartitionWeightedBuffer[ADAPTIVE_PARTITION_MAX_NUM_CANDS]; //Ԥ��ֵ�Ļ��� ������6*5*���ַ�ʽ
	static_vector<uint8_t, ADAPTIVE_PARTITION_MAX_NUM_CANDS> adpPartitionRdModeList; //������Ԥ���RD��ģʽ�б������ڴ�ѡѡ������Ϻ�ѡ���뵽���յ�RDCost��ѡ�б���
	static_vector<double, ADAPTIVE_PARTITION_MAX_NUM_CANDS> adpPartitionCandCostList; //������Ԥ��ĺ�ѡ�����б������ڴ�ѡѡ������Ϻ�ѡ���뵽���յĺ�ѡ�б���

																					  //slice.getMaxNumTriangleCand()����������
																					  // uint8_t numTriangleCandComb = slice.getMaxNumTriangleCand() * (slice.getMaxNumTriangleCand() - 1) * 2;//����Ԥ��ģʽ�����������6*5*2=60 (ֱ�Ӱ�ģ��ø�) 
	uint8_t numTriangleCandComb = 30;

	DistParam distParam;
	const bool useHadamard = !tempCS->slice->getDisableSATDForRD(); //�Ƿ�HAD-Cost�� VTM6.0�ж���!encTestMode.lossless&&!tempCS->slice->getDisableSATDForRD()
	m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, useHadamard);

	//���ٱ���CU�Ŀռ�
	const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));

	const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda();

	MergeCtx mergeCtx; //����AGPM����Ԥ���б����ӳ���Merge�б�����������
	{
		CodingUnit cu(tempCS->area);
		cu.cs = tempCS;
		cu.predMode = MODE_INTER; //֡�����
		cu.slice = tempCS->slice; //ָ��ǰslice
		cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
		cu.adpPartition = true;
		cu.mmvdSkip = false;
		// cu.GBiIdx   = GBI_DEFAULT;//˫��Ԥ��Ȩ������ 6.0��д��
		//�����Ǵ�Geo�︴�Ƶģ��д�����
		cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
		cu.qp = encTestMode.qp;
		cu.affine = false;
		cu.mtsFlag = false;
		cu.bcwIdx = BCW_DEFAULT; //BCW ��˫��Ԥ�⼼��
		cu.geoFlag = false; //GPM��Ϊfalse
		cu.imv = 0;
		cu.mmvdSkip = false;
		cu.skip = false;
		cu.mipFlag = false;
		cu.bdpcmMode = 0;

		PredictionUnit pu(tempCS->area); //Ϊ��ǰ��CU������ʱ����ռ�    
		pu.cu = &cu;
		pu.cs = tempCS;
		pu.regularMergeFlag = false; //����regularģʽ

		//��ȡ��ѡ
		PU::getAGPMMergeCandidates(pu, mergeCtx);

		//const uint8_t maxNumTriangleCand = pu.cs->slice->getMaxNumTriangleCand();//����Ԥ�������ѡ��
		const uint8_t maxNumAdpPartitionCand = 6; //����Ԥ�������ѡ�� ͬ��

												  //����6������Merge��ѡ
		for (uint8_t mergeCand = 0; mergeCand < maxNumAdpPartitionCand; mergeCand++)
		{
			adpPartitionBuffer[mergeCand] = m_acMergeBuffer[mergeCand].getBuf(localUnitArea);
			mergeCtx.setMergeInfo(pu, mergeCand); //����AGPM��ѡ�ĳ�ʼ��Ϣ
			PU::spanMotionInfo(pu, mergeCtx);

			if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(pu))))
			{
				//��ʹ�ô�ģʽ
				tempCS->initStructData(encTestMode.qp); //��ʼ�����ݽṹ
				return;
			}
			m_pcInterSearch->motionCompensation(pu, adpPartitionBuffer[mergeCand]); //�˶����������ÿ����ѡ�ĵ���Ԥ��ֵ
		}
	}

	adpPartitionNumMrgSATDCand = std::min(adpPartitionNumMrgSATDCand, numAdpPartitionCand); //��Ҫ����SADTCost�Ƚϵ�AGPM��ѡ����

	{
		CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType); //����CU

		partitioner.setCUData(cu); //��ʼ��cu
		cu.slice = tempCS->slice;
		cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
		cu.skip = false;
		cu.predMode = MODE_INTER;
		cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
		cu.qp = encTestMode.qp;
		cu.geoFlag = false;
		cu.mmvdSkip = false;
		cu.bcwIdx = BCW_DEFAULT; //˫��Ԥ��Ȩ������ΪĬ��ֵ
		cu.adpPartition = true;

		PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType); //����PU
		numAdpPartitionCand = ADAPTIVE_PARTITION_MAX_NUM_CANDS;
		//numAdpPartitionCand = min(numAdpPartitionCand, numTriangleCandComb);

		//����30�ֺ�ѡ������SAD Cost���д�ѡ
		for (uint8_t mergecand = 0; mergecand < numAdpPartitionCand; mergecand++)
		{
			uint8_t candIdx0 = m_adpPartitionModeTest[mergecand].m_candIdx0; //����һ�ĺ�ѡ����
			uint8_t candIdx1 = m_adpPartitionModeTest[mergecand].m_candIdx1;

			//������ı���������ǰPU
			pu.adpPartitionMergeIdx0 = candIdx0;
			pu.adpPartitionMergeIdx1 = candIdx1;

			pu.mergeFlag = true;
			pu.regularMergeFlag = false;

			adpPartitionWeightedBuffer[mergecand] = m_adpPartitionWeightedBuffer[mergecand].getBuf(localUnitArea); //����Ȩ�ص���ʱ����?
			adpPartitionBuffer[candIdx0] = m_acMergeBuffer[candIdx0].getBuf(localUnitArea); //��һ��������Ԥ��ֵ����ʱ����?
			adpPartitionBuffer[candIdx1] = m_acMergeBuffer[candIdx1].getBuf(localUnitArea);

			//���ȷ�����Ȩ����
			m_pcInterSearch->weightedAdpPartitionBlk(pu, CHANNEL_TYPE_LUMA, adpPartitionWeightedBuffer[mergecand], adpPartitionBuffer[candIdx0], adpPartitionBuffer[candIdx1]);

			distParam.cur = adpPartitionWeightedBuffer[mergecand].Y(); //��ǰ���ģʽ����������Ȩ��
			
			Distortion uiSad = distParam.distFunc(distParam); //���㵱ǰ���ģʽʧ���SAD
			
			uint32_t uiBitsCand = m_adpPartitionIdxBins[candIdx0][candIdx1]; //��ǰ���ģʽ���ܱ��������в�+ͷ��Ϣ+�﷨Ԫ�أ�

			double cost = (double)uiSad + (double)uiBitsCand*sqrtLambdaForFirstPass; //����SAD����

			updateCandList(mergecand, cost, adpPartitionRdModeList, adpPartitionCandCostList, adpPartitionNumMrgSATDCand); //����RD��ѡ�б�����cost��С�ķ����б���
		}
		//����������֮��triangleRdModeList�б��еĺ�ѡ��Cost�Ǵ�С��������ģ�RD�б�����󳤶�Ϊ3��Ҳ���Ǵ�60�����RDCost��3�����ŵĺ�ѡ

		// limit number of candidates using SATD-costs
		//ʹ��SATD-cost���ƺ�ѡ�б�������������һ��SAD��ѡ��֮���RD��ѡ�б����������һ�������ƣ�����RD��ѡ�б��ķ�Χ����3��������������������Ϊ0��1��2��3������һ������
		for (uint8_t i = 0; i < adpPartitionNumMrgSATDCand; i++)
		{
			//�����ǰRD��ѡ��Cost����һ������ֵ�������������ǰ����Щ��Cost��ֵ��Χ֮�ڵĺ�ѡ�������б��У�������Щ����Cost��ֵ�ĺ�ѡ���ų�����
			if (adpPartitionCandCostList[i] > MRG_FAST_RATIO * adpPartitionCandCostList[0] || adpPartitionCandCostList[i] > getMergeBestSATDCost())
			{
				adpPartitionNumMrgSATDCand = i;//STAD��ѡ����������Ϊ��Щ����ֵ֮�ڵĺ�ѡ����
				break;
			}
		}

		// perform chroma weighting process
		//ִ��ɫ�ȼ�Ȩ����

		//����Щ����SAD�Լ���һ����������ɸѡ���Ľ��ŵĺ�ѡ���б�����ע�⣺�����ȷ�������ɸѡ��ɸѡ��֮����ȥ����Ӧ��ɫ�Ƚ��м�Ȩ��������Ϊ���ǵ��������ԭ��һ�����ȷ����Ŀɿ��Ը�ǿ��������ɸѡ���ĺ�ѡ�Ƚ���ʵ�ɿ�������Ϊ�˽�ʡ���Ӷȣ�ֻ��Ҫ�����ȷ���ȥ����ɸѡ���㹻�ˣ�����Ҫ��ȥ��ɫ�ȷ�����
		for (uint8_t i = 0; i < adpPartitionNumMrgSATDCand; i++)
		{
			uint8_t  mergeCand = adpPartitionRdModeList[i]; //RdModeList�б��еĵڼ�����ѡ
			uint8_t  candIdx0 = m_adpPartitionModeTest[mergeCand].m_candIdx0;//����һ�ĺ�ѡ
			uint8_t  candIdx1 = m_adpPartitionModeTest[mergeCand].m_candIdx1;//�������ĺ�ѡ

			pu.adpPartitionMergeIdx0 = candIdx0;
			pu.adpPartitionMergeIdx1 = candIdx1;
			pu.mergeFlag = true;//��־Mergeģʽ����
			pu.regularMergeFlag = false;//��־regular_Mergeģʽ������

										//ɫ�ȷ�����Ȩ����
			m_pcInterSearch->weightedAdpPartitionBlk(pu, CHANNEL_TYPE_CHROMA, adpPartitionWeightedBuffer[mergeCand], adpPartitionBuffer[candIdx0], adpPartitionBuffer[candIdx1]);
		}

		//�ٳ�ʼ�����ݽṹ
		tempCS->initStructData(encTestMode.qp);
	}
	//���»�ȡ����Ԥ��ģʽ������������֮���triangleNumMrgSATDCand��RD�б��ĳߴ���ȡ��С
	adpPartitionNumMrgSATDCand = std::min(adpPartitionNumMrgSATDCand, (uint8_t)adpPartitionRdModeList.size());

	m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;//�Ƿ������ŵ�Cost


																	  //���￪ʼ��ʽ��HAD_Cost��ϸѡ����
	{
		uint8_t iteration;
		uint8_t iterationBegin = 0;//������0��ʼ
		/*
		if (encTestMode.lossless)//������Եĺ�ѡģʽ���������һ��
		{
			iteration = 1;
		}
		else//�����������
		{
			iteration = 2;
		}
		*/
		iteration = 2;

		//���ĵ�������
		for (uint8_t noResidualPass = iterationBegin; noResidualPass < iteration; ++noResidualPass)
		{
			//������Ҫ����SATDϸѡ���б��е�ÿһ����ѡ������SATD_Cost�ıȽϣ��õ�cost��С���Ǹ���ѡ��Ϊ���ŵ�����ģʽ��Ȼ����ȥ�������Mergeģʽȥ������ʧ����۵ıȽ�
			for (uint8_t mrgHADIdx = 0; mrgHADIdx < adpPartitionNumMrgSATDCand; mrgHADIdx++)
			{
				uint8_t mergeCand = adpPartitionRdModeList[mrgHADIdx];

				if (((noResidualPass != 0) && AdpPartitionHasNoResidual[mergeCand])
					|| ((noResidualPass == 0) && bestIsSkip))
				{
					continue;
				}

				uint8_t candIdx0 = m_adpPartitionModeTest[mergeCand].m_candIdx0;
				uint8_t candIdx1 = m_adpPartitionModeTest[mergeCand].m_candIdx1;

				CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);

				partitioner.setCUData(cu);//����CU������
				cu.slice = tempCS->slice;
				cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
				cu.skip = false;
				cu.predMode = MODE_INTER;
				cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
				cu.qp = encTestMode.qp;
				cu.geoFlag = false;
				cu.mmvdSkip = false;
				cu.bcwIdx = BCW_DEFAULT; //˫��Ԥ��Ȩ������ΪĬ��ֵ
				cu.adpPartition = true;
				cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
				cu.qp = encTestMode.qp;
				cu.mmvdSkip = false;
				PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);

				pu.adpPartitionMergeIdx0 = candIdx0;
				pu.adpPartitionMergeIdx1 = candIdx1;
				pu.mergeFlag = true;
				pu.regularMergeFlag = false;

				//������������Ԥ���˶���Ϣ
				PU::spanAGPMMotionInfo(pu, mergeCtx, candIdx0, candIdx1);

				if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(*cu.firstPU))))
				{
					// Do not use this mode
					tempCS->initStructData(encTestMode.qp);
					return;
				}

				//��֮ǰ�õ��ļ�Ȩ���Ԥ��ֱֵ�ӿ�����������������Ĳв���루���оͽ�����HAD��cost�ļ����Լ��Ƚϣ�
				tempCS->getPredBuf().copyFrom(adpPartitionWeightedBuffer[mergeCand]);//�õ�����Ԥ�������Ԥ��ֵ

																				 //����֡��в�
				xEncodeInterResidual(tempCS, bestCS, partitioner, encTestMode, noResidualPass, (noResidualPass == 0 ? &AdpPartitionHasNoResidual[mergeCand] : NULL));

				if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
				{
					bestIsSkip = bestCS->getCU(partitioner.chType)->rootCbf == 0;
				}
				tempCS->initStructData(encTestMode.qp);
			}// end loop mrgHADIdx
		}
	}
	if (m_bestModeUpdated && bestCS->cost != MAX_DOUBLE)//ѡ�����ŵ�һ���������ģʽ���������������Mergeģʽһ��ȥ����
	{
		xCalDebCost(*bestCS, partitioner);
	}
}

void PU::spanAGPMMotionInfo(PredictionUnit &pu, MergeCtx &mrgCtx, const uint8_t candIdx0, const uint8_t candIdx1)
{
	pu.adpPartitionMergeIdx0 = candIdx0;
	pu.adpPartitionMergeIdx1 = candIdx1;
	MotionBuf mb = pu.getMotionBuf();

	MotionInfo biMv;
	biMv.isInter = true;
	biMv.sliceIdx = pu.cs->slice->getIndependentSliceIdx();

	if (mrgCtx.interDirNeighbours[candIdx0] == 1 && mrgCtx.interDirNeighbours[candIdx1] == 2)
	{
		biMv.interDir = 3;
		biMv.mv[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].mv;
		biMv.mv[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
		biMv.refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx;
		biMv.refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
	}
	else if (mrgCtx.interDirNeighbours[candIdx0] == 2 && mrgCtx.interDirNeighbours[candIdx1] == 1)
	{
		biMv.interDir = 3;
		biMv.mv[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
		biMv.mv[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
		biMv.refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx;
		biMv.refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
	}
	else if (mrgCtx.interDirNeighbours[candIdx0] == 1 && mrgCtx.interDirNeighbours[candIdx1] == 1)
	{
#if AGPM_SIMPLIFICATION	
		biMv.interDir = 1;
		biMv.mv[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
		biMv.mv[1] = Mv(0, 0);
		biMv.refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx;
		biMv.refIdx[1] = -1;
#else
		int32_t refIdx = mappingRefPic(pu, pu.cs->slice->getRefPOC(REF_PIC_LIST_0, mrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx), REF_PIC_LIST_1);
		if (refIdx != -1)
		{
			biMv.interDir = 3;
			biMv.mv[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].mv;
			biMv.mv[1] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
			biMv.refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx;
			biMv.refIdx[1] = refIdx;
		}
		else
		{
			refIdx = mappingRefPic(pu, pu.cs->slice->getRefPOC(REF_PIC_LIST_0, mrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx), REF_PIC_LIST_1);
			biMv.interDir = (refIdx != -1) ? 3 : 1;
			biMv.mv[0] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[candIdx1 << 1].mv : mrgCtx.mvFieldNeighbours[candIdx0 << 1].mv;
			biMv.mv[1] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[candIdx0 << 1].mv : Mv(0, 0);
			biMv.refIdx[0] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx : mrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx;
			biMv.refIdx[1] = (refIdx != -1) ? refIdx : -1;
		}
#endif
	}
	else if (mrgCtx.interDirNeighbours[candIdx0] == 2 && mrgCtx.interDirNeighbours[candIdx1] == 2)
	{
#if AGPM_SIMPLIFICATION	
		biMv.interDir = 2;
		biMv.mv[0] = Mv(0, 0);
		biMv.mv[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
		biMv.refIdx[0] = -1;
		biMv.refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
#else
		int32_t refIdx = mappingRefPic(pu, pu.cs->slice->getRefPOC(REF_PIC_LIST_1, mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx), REF_PIC_LIST_0);
		if (refIdx != -1)
		{
			biMv.interDir = 3;
			biMv.mv[0] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
			biMv.mv[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
			biMv.refIdx[0] = refIdx;
			biMv.refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
		}
		else
		{
			refIdx = mappingRefPic(pu, pu.cs->slice->getRefPOC(REF_PIC_LIST_1, mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx), REF_PIC_LIST_0);
			biMv.interDir = (refIdx != -1) ? 3 : 2;
			biMv.mv[0] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv : Mv(0, 0);
			biMv.mv[1] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv : mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
			biMv.refIdx[0] = (refIdx != -1) ? refIdx : -1;
			biMv.refIdx[1] = (refIdx != -1) ? mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx : mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
		}
#endif
	}

	int32_t idxW = (int32_t)(floorLog2(pu.lwidth()) - MIN_CU_LOG2);
	int32_t idxH = (int32_t)(floorLog2(pu.lheight()) - MIN_CU_LOG2);
	for (int32_t y = 0; y < mb.height; y++)
	{
		for (int32_t x = 0; x < mb.width; x++)
		{
			if (g_adpPartitionMvStorage[idxH][idxW][y][x] == 2)
			{
				mb.at(x, y).isInter = true;
				mb.at(x, y).interDir = biMv.interDir;
				mb.at(x, y).refIdx[0] = biMv.refIdx[0];
				mb.at(x, y).refIdx[1] = biMv.refIdx[1];
				mb.at(x, y).mv[0] = biMv.mv[0];
				mb.at(x, y).mv[1] = biMv.mv[1];
				mb.at(x, y).sliceIdx = biMv.sliceIdx;
			}
			else if (g_adpPartitionMvStorage[idxH][idxW][y][x] == 0)
			{
				mb.at(x, y).isInter = true;
				mb.at(x, y).interDir = mrgCtx.interDirNeighbours[candIdx0];
				mb.at(x, y).refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].refIdx;
				mb.at(x, y).refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].refIdx;
				mb.at(x, y).mv[0] = mrgCtx.mvFieldNeighbours[candIdx0 << 1].mv;
				mb.at(x, y).mv[1] = mrgCtx.mvFieldNeighbours[(candIdx0 << 1) + 1].mv;
				mb.at(x, y).sliceIdx = biMv.sliceIdx;
			}
			else
			{
				mb.at(x, y).isInter = true;
				mb.at(x, y).interDir = mrgCtx.interDirNeighbours[candIdx1];
				mb.at(x, y).refIdx[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].refIdx;
				mb.at(x, y).refIdx[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].refIdx;
				mb.at(x, y).mv[0] = mrgCtx.mvFieldNeighbours[candIdx1 << 1].mv;
				mb.at(x, y).mv[1] = mrgCtx.mvFieldNeighbours[(candIdx1 << 1) + 1].mv;
				mb.at(x, y).sliceIdx = biMv.sliceIdx;
			}
		}
	}
}   

#if !AGPM_SIMPLIFICATION
int32_t PU::mappingRefPic(const PredictionUnit &pu, int32_t refPicPoc, bool targetRefPicList)
{
	int32_t numRefIdx = pu.cs->slice->getNumRefIdx((RefPicList)targetRefPicList);

	for (int32_t i = 0; i < numRefIdx; i++)
	{
		if (pu.cs->slice->getRefPOC((RefPicList)targetRefPicList, i) == refPicPoc)
		{
			return i;
		}
	}
	return -1;
}
#endif

/*
bool PU::isUniqueAdpPartitionCandidates(const PredictionUnit &pu, MergeCtx& mrgCtx)
{
	int newCand = mrgCtx.numValidMergeCand;
	for (int32_t i = 0; i < newCand; i++)
	{
		int32_t predFlagCur = mrgCtx.interDirNeighbours[i] == 1 ? 0 : 1;
		int32_t predFlagNew = mrgCtx.interDirNeighbours[newCand] == 1 ? 0 : 1;
		int32_t refPicPocCur = pu.cs->slice->getRefPOC((RefPicList)predFlagCur, mrgCtx.mvFieldNeighbours[(i << 1) + predFlagCur].refIdx);
		int32_t refPicPocNew = pu.cs->slice->getRefPOC((RefPicList)predFlagNew, mrgCtx.mvFieldNeighbours[(newCand << 1) + predFlagNew].refIdx);
		if (refPicPocCur == refPicPocNew && mrgCtx.mvFieldNeighbours[(i << 1) + predFlagCur].mv == mrgCtx.mvFieldNeighbours[(newCand << 1) + predFlagNew].mv)
		{
			return false;
		}
	}
	return true;
}
*/
#endif