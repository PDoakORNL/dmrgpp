/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file DensityMatrixSu2.h
 *
 *
 *
 */
#ifndef DENSITY_MATRIX_SU2_H
#define DENSITY_MATRIX_SU2_H

#include "Matrix.h" // in PsimagLite
#include "AlmostEqual.h" // in PsimagLite
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"

namespace Dmrg {
	template<typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
	class DensityMatrixSu2 : public DensityMatrixBase<DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType> {
		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef BlockMatrix<PsimagLite::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		typedef typename PsimagLite::Real<DensityMatrixElementType>::Type RealType;

		enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

	public:
		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;

		DensityMatrixSu2(
			const TargettingType&,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType&,
			const DmrgBasisType&,
			SizeType,bool debug=false,bool verbose=false) : data_(pBasis.size() ,pBasis.partition()-1),mMaximal_(pBasis.partition()-1),pBasis_(pBasis),
				debug_(debug),verbose_(verbose)
		{
		}

		BlockMatrixType& operator()()
		{
			return data_;
		}

		SizeType rank() { return data_.rank(); }

		void check(int direction)
		{
			if (!debug_) return;

			if (verbose_) std::cerr<<"CHECKING DENSITY-MATRIX WITH OPTION="<<direction<<"\n";
			for (SizeType m=0;m<data_.blocks();m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
				SizeType p = mMaximal_[m];
				if (m==p) continue;
				//is data_(m)==data_(p) ?
				check(m,data_(m),p,data_(p));

			}
		}

		void check2(int direction)
		{
			if (!debug_) return;
			if (verbose_) std::cerr<<"CHECKING DMRG-TRANFORM WITH OPTION="<<direction<<"\n";
			for (SizeType m=0;m<data_.blocks();m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
				SizeType p = mMaximal_[m];
				if (m==p) continue;
				//is data_(m)==data_(p) ?
				check2(m,data_(m),p,data_(p));

			}
		}

		void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
		{
			diagonalise(data_,eigs,jobz);

			//make sure non-maximals are equal to maximals
			// this is needed because otherwise there's no assure that m-independence
			// is achieved due to the non unique phase of eigenvectors of the density matrix
			for (SizeType m=0;m<data_.blocks();m++) {

				SizeType p = mMaximal_[m];
				if (SizeType(m)==p) continue; // we already did these ones

				data_.setBlock(m,data_.offsets(m),data_(p));
			}
			if (verbose_) std::cerr<<"After diagonalise\n";

			if (debug_) areAllMsEqual(pBasis_);

		}

		void init(
				const TargettingType& target,
				DmrgBasisWithOperatorsType const &pBasis,
				const DmrgBasisWithOperatorsType& pBasisSummed,
				DmrgBasisType const &pSE,
				int direction)
		{
			BuildingBlockType matrixBlock;

			for (SizeType m=0;m<pBasis.partition()-1;m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)

				if (DmrgBasisType::useSu2Symmetry()) {
					mMaximal_[m] = findMaximalPartition(m,pBasis);
					//if (enforceSymmetry && SizeType(m)!=mMaximal_[m]) continue;
					// we'll fill non-maximal partitions later
				}

				SizeType bs = pBasis.partition(m+1)-pBasis.partition(m);

				matrixBlock.reset(bs,bs);

				for (SizeType i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
					for (SizeType j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {

						matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m))=
							densityMatrixAux(i,j,target,pBasisSummed,pSE,direction);

					}
				}
				data_.setBlock(m,pBasis.partition(m),matrixBlock);
			}

			if (verbose_) {
				std::cerr<<"DENSITYMATRIXPRINT option="<<direction<<"\n";
				std::cerr<<(*this);
				std::cerr<<"***********\n";
				std::cerr<<"Calling ae from init()...\n";
			}
			if (debug_) areAllMsEqual(pBasis);
		}

		template<typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			>
		friend std::ostream& operator<<(std::ostream& os,
				const DensityMatrixSu2<
    					DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_>& dm);
	private:
		BlockMatrixType data_;
		typename PsimagLite::Vector<SizeType>::Type mMaximal_;
		const DmrgBasisWithOperatorsType& pBasis_;
		bool debug_,verbose_;

		SizeType findMaximalPartition(SizeType p,DmrgBasisWithOperatorsType const &pBasis)
		{
			std::pair<SizeType,SizeType> jm2 = pBasis.jmValue(pBasis.partition(p));
			SizeType ne2 = pBasis.electrons(pBasis.partition(p));
			if (jm2.first==jm2.second) return p;
			for (SizeType m=0;m<pBasis.partition()-1;m++) {
				std::pair<SizeType,SizeType> jm1 = pBasis.jmValue(pBasis.partition(m));
				SizeType ne1 = pBasis.electrons(pBasis.partition(m));
				if (jm1.first==jm2.first && jm1.first==jm1.second && ne1==ne2) return m;
			}
			throw PsimagLite::RuntimeError("	findMaximalPartition : none found\n");
		}

		//! only used for debugging
		bool areAllMsEqual(DmrgBasisWithOperatorsType const &)
		{
//			PsimagLite::AlmostEqual<RealType> almostEqual(1e-5);

//			for (SizeType m=0;m<pBasis.partition()-1;m++) {
//				std::pair<SizeType,SizeType> jmPair1 =  pBasis.jmValue(pBasis.partition(m));
//				SizeType ne1 = pBasis.electrons(pBasis.partition(m));
//				for (SizeType p=0;p<pBasis.partition()-1;p++) {
//					std::pair<SizeType,SizeType> jmPair2 =  pBasis.jmValue(pBasis.partition(p));
//					SizeType ne2 = pBasis.electrons(pBasis.partition(p));

//					if (jmPair1.first == jmPair2.first && ne1==ne2) {

//						if (!almostEqual(data_(m),data_(p))) {
//							std::cerr<<"Checking m="<<m<<" against p="<<p<<"\n";
//							throw PsimagLite::RuntimeError("error\n");
//						}
//					}
//				}
//			}
			return true;
		}
		DensityMatrixElementType densityMatrixAux(SizeType alpha1,SizeType alpha2,const TargettingType& target,
			DmrgBasisWithOperatorsType const &pBasisSummed,DmrgBasisType const &pSE,SizeType direction)
		{
			DensityMatrixElementType sum=0;
			// The g.s. has to be treated separately because it's usually a vector of RealType, whereas
			// the other targets might be complex, and C++ generic programming capabilities are weak... we need D!!!
			if (target.includeGroundStage()) sum +=  densityMatrixHasFactors(alpha1,alpha2,target.gs(),
			    		pBasisSummed,pSE,direction)*target.gsWeight();

			for (SizeType i=0;i<target.size();i++)
				sum += densityMatrixHasFactors(alpha1,alpha2,target(i),
					pBasisSummed,pSE,direction)*target.weight(i)/target.normSquared(i);

			return sum;
		}

		template<typename TargetVectorType>
		DensityMatrixElementType densityMatrixHasFactors(SizeType alpha1,SizeType alpha2,const TargetVectorType& v,
			DmrgBasisWithOperatorsType const &pBasisSummed,DmrgBasisType const &pSE,SizeType direction)
		{
			int ne = pBasisSummed.size();
			int ns = pSE.size()/ne;
			SizeType total=pBasisSummed.size();
			if (direction!=EXPAND_SYSTEM) {
				ns=pBasisSummed.size();
				ne=pSE.size()/ns;
			}

			DensityMatrixElementType sum=0;

			// Make sure we don't copy just get the reference here!!
			const FactorsType& factors = pSE.getFactors();

			for (SizeType beta=0;beta<total;beta++) {
				// sum over environ:
				int i1 = alpha1+beta*ns;
				int i2 = alpha2+beta*ns;
				// sum over system:
				if (direction!=EXPAND_SYSTEM) {
					i1 = beta + alpha1*ns;
					i2 = beta + alpha2*ns;
				}
				for (int k1=factors.getRowPtr(i1);k1<factors.getRowPtr(i1+1);k1++) {
					int eta1 = factors.getCol(k1);
					int ii = pSE.permutationInverse(eta1);
					for (int k2=factors.getRowPtr(i2);k2<factors.getRowPtr(i2+1);k2++) {
						int eta2 =  factors.getCol(k2);
						int jj = pSE.permutationInverse(eta2);

						DensityMatrixElementType tmp3= v.slowAccess(ii)*
						        PsimagLite::conj(v.slowAccess(jj)) *
							factors.getValue(k1) * factors.getValue(k2);
						sum += tmp3;
					}
				}
			}
			return sum;
		}

		//! only used for debugging
		void check(SizeType p1,const BuildingBlockType& bp1,SizeType p2,const BuildingBlockType& bp2)
		{
			if (bp1.n_row()!=bp2.n_row()) {
				std::cerr<<"row size different "<<bp1.n_row()<<"!="<<bp2.n_row()<<"\n";
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				std::cerr<<data_<<"\n";
				throw PsimagLite::RuntimeError("Density Matrix Check: failed\n");
			}
			if (bp1.n_col()!=bp2.n_col()) {
				std::cerr<<"col size different "<<bp1.n_col()<<"!="<<bp2.n_col()<<"\n";
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				std::cerr<<data_<<"\n";
				throw PsimagLite::RuntimeError("Density Matrix Check: failed\n");
			}
			if (!debug_) return;

			for (SizeType i=0;i<bp1.n_row();i++) {
				for (SizeType j=0;j<bp1.n_col();j++) {
					RealType x = PsimagLite::norm(bp1(i,j)-bp2(i,j));
					if (x>1e-5) {
						std::cerr<<bp1(i,j)<<"!="<<bp2(i,j)<<" i="<<i<<" j= "<<j<<"\n";
						std::cerr<<"difference="<<x<<" p1="<<p1<<" p2="<<p2<<"\n";
						std::cerr<<data_(p1)<<"\n";
						std::cerr<<"******************\n";
						std::cerr<<data_(p2)<<"\n";
						throw PsimagLite::RuntimeError("Density Matrix Check: failed (differ)\n");
					}
				}
			}
		}

		//! only used for debugging
		void check2(SizeType p1,const BuildingBlockType& bp1,SizeType p2,const BuildingBlockType& bp2)
		{
			DensityMatrixElementType alpha=1.0,beta=0.0;
			int n =bp1.n_col();
			PsimagLite::Matrix<DensityMatrixElementType> result(n,n);
			psimag::BLAS::GEMM('C','N',n,n,n,alpha,&(bp1(0,0)),n,&(bp2(0,0)),n,beta,&(result(0,0)),n);
			if (!PsimagLite::isTheIdentity(result)) {
				//utils::matrixPrint(result,std::cerr);
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				throw PsimagLite::RuntimeError("Density Matrix Check2: failed\n");
			}

		}

	}; // class DensityMatrixSu2

	template<typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		>
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrixSu2<DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		//std::cerr<<"PRINTING DENSITY-MATRIX WITH OPTION="<<option<<"\n";
		for (SizeType m=0;m<dm.data_.blocks();m++) {
			// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
			std::pair<SizeType,SizeType> jm1 = dm.pBasis_.jmValue(dm.pBasis_.partition(m));
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<"partitionNumber="<<m<<" j="<<jm1.first<<" m= "<<jm1.second<<" ne="<<ne<<"\n";
			os<<dm.data_(m)<<"\n";
		}
		return os;
	}
} // namespace Dmrg

/*@}*/
#endif


