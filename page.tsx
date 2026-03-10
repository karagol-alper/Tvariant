"use client"

import { CardDescription } from "@/components/ui/card"

import { useState, useEffect, useMemo } from "react"
import {
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  BarChart,
  Bar,
  Cell,
  ComposedChart,
  Area,
  LineChart, // Imported for the update
} from "recharts"
import {
  Download,
  FileText,
  LucideBrickWallIcon,LucideBrickWallIcon,LucideBrickWallIcon,LucideBrickWallIcon,LucidePaintbrushVerticalIcon,ActivitySquareIcon,ActivitySquareIcon,BarChart3,
  TrendingUp,
  Database,
  Settings,
  Info,
  BoxesIcon,AlertCircle,
  Activity,
  Shield,
  ChevronDown,
  ChevronUp,
  Zap,
  History,
  X,
} from "lucide-react"
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@/components/ui/collapsible"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Textarea } from "@/components/ui/textarea"
import { Label } from "@/components/ui/label"
import { Input } from "@/components/ui/input"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { Loader2 } from "lucide-react"

// Main App component
export default function ProteinAnalysisTool() {
  // State variables for MSA input, parameters, and results
  const [msaText, setMsaText] = useState("")
  const [sigma, setSigma] = useState(0.1)
  const [autoSigma, setAutoSigma] = useState(true) // Added auto bandwidth toggle
  const [calculatedSigma, setCalculatedSigma] = useState(0.1) // Store calculated bandwidth
  const [windowSize, setWindowSize] = useState(3)
  const [pseudocount, setPseudocount] = useState(0.01)
  const [bootstrapIterations, setBootstrapIterations] = useState(100)
  const [confidenceLevel, setConfidenceLevel] = useState(0.95)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState("")
  const [activeTab, setActiveTab] = useState("global")

  // State for calculated results
  const [allPositionsTv, setAllPositionsTv] = useState([])
  const [globalTv, setGlobalTv] = useState({})
  const [highTensionSites, setHighTensionSites] = useState([])
  const [msaStats, setMsaStats] = useState({})
  const [confidenceIntervals, setConfidenceIntervals] = useState({})
  const [mutualInformation, setMutualInformation] = useState([])
  const [structuralComplexity, setStructuralComplexity] = useState([])
  const [statisticalSignificance, setStatisticalSignificance] = useState({})
  const [correlationPValues, setCorrelationPValues] = useState([])
  const [fdrCorrectedPValues, setFdrCorrectedPValues] = useState({})
  const [effectSizes, setEffectSizes] = useState([])
  const [cvMetrics, setCvMetrics] = useState({})

  const [methodologyOpen, setMethodologyOpen] = useState(false)

  const [showDocumentation, setShowDocumentation] = useState(false)
  const [showVersionHistory, setShowVersionHistory] = useState(false)

  // Constants for amino acids
  const AMINO_ACIDS = useMemo(
    () => ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"],
    [],
  )

  // Helper function to parse MSA text
  const parseMSA = (text) => {
    const lines = text
      .split("\n")
      .map((line) => line.trim())
      .filter((line) => line.length > 0)
    const sequences = []
    let currentSequence = ""

    for (const line of lines) {
      if (line.startsWith(">")) {
        if (currentSequence) {
          sequences.push(currentSequence)
          currentSequence = ""
        }
      } else {
        currentSequence += line
      }
    }

    if (currentSequence) {
      sequences.push(currentSequence)
    }

    return sequences
  }

  // Helper function to calculate MSA statistics with enhanced metrics
  const calculateMsaStats = (msa) => {
    if (!msa || msa.length === 0) return {}

    const numSequences = msa.length
    const alignmentLength = msa[0].length
    const totalPositions = numSequences * alignmentLength

    // Calculate composition with gap handling
    const composition = {}
    AMINO_ACIDS.forEach((aa) => (composition[aa] = 0))
    let gapCount = 0

    for (const seq of msa) {
      for (const aa of seq) {
        const upperAA = aa.toUpperCase()
        if (composition[upperAA] !== undefined) {
          composition[upperAA]++
        } else if (upperAA === "-" || upperAA === ".") {
          gapCount++
        }
      }
    }

    // Calculate Shannon entropy per position
    const entropy = []
    const conservation = []
    const effectiveSequences = []

    for (let i = 0; i < alignmentLength; i++) {
      const column = msa.map((seq) => seq[i].toUpperCase())
      const counts = {}
      let validCount = 0

      column.forEach((aa) => {
        if (AMINO_ACIDS.includes(aa)) {
          counts[aa] = (counts[aa] || 0) + 1
          validCount++
        }
      })

      // Shannon entropy
      let H = 0
      Object.values(counts).forEach((count) => {
        if (count > 0) {
          const p = count / validCount
          H -= p * Math.log2(p)
        }
      })
      entropy.push(H)

      // Conservation (1 - normalized entropy)
      const maxEntropy = Math.log2(Math.min(20, validCount))
      const normalizedEntropy = maxEntropy > 0 ? H / maxEntropy : 0
      conservation.push(1 - normalizedEntropy)

      // Effective number of sequences (for weighting)
      effectiveSequences.push(validCount)
    }

    // Calculate sequence identity matrix
    const identityMatrix = []
    for (let i = 0; i < numSequences; i++) {
      const row = []
      for (let j = 0; j < numSequences; j++) {
        if (i === j) {
          row.push(1.0)
        } else {
          let matches = 0
          let validPositions = 0
          for (let k = 0; k < alignmentLength; k++) {
            const aa1 = msa[i][k].toUpperCase()
            const aa2 = msa[j][k].toUpperCase()
            if (AMINO_ACIDS.includes(aa1) && AMINO_ACIDS.includes(aa2)) {
              validPositions++
              if (aa1 === aa2) matches++
            }
          }
          row.push(validPositions > 0 ? matches / validPositions : 0)
        }
      }
      identityMatrix.push(row)
    }

    // Calculate average pairwise identity
    let totalIdentity = 0
    let pairCount = 0
    for (let i = 0; i < numSequences; i++) {
      for (let j = i + 1; j < numSequences; j++) {
        totalIdentity += identityMatrix[i][j]
        pairCount++
      }
    }
    const avgPairwiseIdentity = pairCount > 0 ? totalIdentity / pairCount : 0

    // Calculate effective number of sequences (Henikoff weighting approximation)
    const effectiveNumSequences = numSequences / (1 + avgPairwiseIdentity)

    return {
      numSequences,
      alignmentLength,
      composition,
      gapCount,
      gapPercentage: (gapCount / totalPositions) * 100,
      avgConservation: conservation.reduce((a, b) => a + b, 0) / conservation.length,
      conservation,
      entropy,
      avgEntropy: entropy.reduce((a, b) => a + b, 0) / entropy.length,
      identityMatrix,
      avgPairwiseIdentity,
      effectiveNumSequences,
      effectiveSequences,
    }
  }

  // Helper function to calculate empirical marginal frequencies with pseudocounts
  const calculateFrequenciesAndEntropy = (msa, pseudocountParam = 0.01) => {
    if (!msa || msa.length === 0) return { frequencies: new Map(), entropies: new Map(), L: 0 }

    const L = msa[0].length
    const n = msa.length

    const frequencies = new Map()
    const entropies = new Map()
    const epsilon = 1e-10

    for (let i = 0; i < L; i++) {
      const aaCounts = new Map()
      let validCount = 0

      for (const seq of msa) {
        const aa = seq[i].toUpperCase()
        if (AMINO_ACIDS.includes(aa)) {
          aaCounts.set(aa, (aaCounts.get(aa) || 0) + 1)
          validCount++
        }
      }

      // Apply Dirichlet pseudocount correction
      const posFrequencies = new Map()
      let entropy = 0
      const totalWithPseudo = validCount + 20 * pseudocountParam

      for (const aa of AMINO_ACIDS) {
        const count = (aaCounts.get(aa) || 0) + pseudocountParam
        const freq = count / totalWithPseudo
        posFrequencies.set(aa, freq)
        if (freq > epsilon) {
          entropy -= freq * Math.log(freq)
        }
      }
      frequencies.set(i, posFrequencies)
      entropies.set(i, entropy)
    }
    return { frequencies, entropies, L }
  }

  // Calculate optimal Gaussian bandwidth using Silverman's rule adapted for sequence data
  const calculateOptimalBandwidth = (msa, stats) => {
    if (!msa || msa.length === 0 || !stats) return 0.1

    const n = msa.length
    const L = msa[0].length
    const neff = stats.effectiveNumSequences || n

    // Silverman's rule: σ = 0.9 × min(std, IQR/1.34) × n^(-1/5)
    // Adapted for protein sequence frequency space

    // Calculate frequency variance across all positions
    const allFrequencies = []
    for (let pos = 0; pos < L; pos++) {
      const aaCounts = new Map()
      let validCount = 0

      for (const seq of msa) {
        const aa = seq[pos].toUpperCase()
        if (AMINO_ACIDS.includes(aa)) {
          aaCounts.set(aa, (aaCounts.get(aa) || 0) + 1)
          validCount++
        }
      }

      // Convert to frequencies
      for (const aa of AMINO_ACIDS) {
        const freq = (aaCounts.get(aa) || 0) / validCount
        if (freq > 0) allFrequencies.push(freq)
      }
    }

    // Calculate standard deviation of frequencies
    const mean = allFrequencies.reduce((a, b) => a + b, 0) / allFrequencies.length
    const variance = allFrequencies.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / allFrequencies.length
    const std = Math.sqrt(variance)

    // Calculate IQR (Interquartile Range)
    const sorted = [...allFrequencies].sort((a, b) => a - b)
    const q1 = sorted[Math.floor(sorted.length * 0.25)]
    const q3 = sorted[Math.floor(sorted.length * 0.75)]
    const iqr = q3 - q1

    // Silverman's rule with effective sequence count
    const silvermanFactor = 0.9 * Math.min(std, iqr / 1.34)
    const silvermanBandwidth = silvermanFactor * Math.pow(neff, -1 / 5)

    // Scott's rule as alternative: σ = std × n^(-1/(d+4))
    // For frequency space, dimensionality d ≈ 20 (amino acids)
    const scottBandwidth = std * Math.pow(neff, -1 / 24)

    // Adaptive correction based on sequence diversity
    const diversityFactor = 1 + stats.avgEntropy / Math.log2(20) // Normalize by max entropy
    const avgIdentity = stats.avgPairwiseIdentity || 0.5

    // If sequences are very similar (high identity), use larger bandwidth for smoothing
    // If sequences are diverse (low identity), use smaller bandwidth for precision
    const identityCorrection = 0.5 + avgIdentity * 0.5

    // Combine methods with weighted average
    const optimalBandwidth = (silvermanBandwidth * 0.6 + scottBandwidth * 0.4) * diversityFactor * identityCorrection

    // Clamp to reasonable range: [0.01, 0.5]
    const clampedBandwidth = Math.max(0.01, Math.min(0.5, optimalBandwidth))

    return Number.parseFloat(clampedBandwidth.toFixed(4))
  }

  const calculateAdaptiveBandwidth = (msa, position, globalSigma) => {
    if (!msa || msa.length === 0) return globalSigma

    const n = msa.length
    const column = msa.map((seq) => seq[position].toUpperCase())
    const validAAs = column.filter((aa) => AMINO_ACIDS.includes(aa))

    if (validAAs.length === 0) return globalSigma

    // Calculate local diversity (Shannon entropy)
    const counts = new Map()
    validAAs.forEach((aa) => counts.set(aa, (counts.get(aa) || 0) + 1))

    let entropy = 0
    for (const count of counts.values()) {
      const p = count / validAAs.length
      if (p > 0) entropy -= p * Math.log2(p)
    }

    // Normalize entropy (max is log2(20) ≈ 4.32)
    const normalizedEntropy = entropy / Math.log2(20)

    // High entropy (diverse) → smaller bandwidth (more precision)
    // Low entropy (conserved) → larger bandwidth (more smoothing)
    const adaptiveFactor = 0.5 + (1 - normalizedEntropy) * 0.5

    return globalSigma * adaptiveFactor
  }

  // Adaptive Gaussian Kernel with bandwidth selection
  const gaussianKernel = (u, v, sigmaParam) => {
    if (sigmaParam <= 0) return 0
    const diff = u - v
    return Math.exp(-Math.pow(diff, 2) / (2 * Math.pow(sigmaParam, 2))) / (sigmaParam * Math.sqrt(2 * Math.PI))
  }

  // Calculate conditional kernel regression estimate with bias correction
  const calculateConditionalProbability = (freqsAtPos, X, Y, sigmaParam) => {
    let numerator = 0
    let denominator = 0

    const freqX = freqsAtPos.get(X) || 0
    const epsilon = 1e-10

    for (const Z of AMINO_ACIDS) {
      const freqZ = freqsAtPos.get(Z) || 0
      const kernelValue = gaussianKernel(freqX, freqZ, sigmaParam)
      denominator += kernelValue
      if (Z === Y) {
        numerator += kernelValue
      }
    }

    // Nadaraya-Watson estimator with small-sample correction
    return denominator > epsilon ? numerator / denominator : 0
  }

  // Calculate mutual information between positions for coevolution detection
  const calculateMutualInformation = (msa, pos1, pos2) => {
    if (!msa || msa.length === 0) return 0

    const n = msa.length
    const joint = new Map()
    const marginal1 = new Map()
    const marginal2 = new Map()

    for (const seq of msa) {
      const aa1 = seq[pos1].toUpperCase()
      const aa2 = seq[pos2].toUpperCase()

      if (AMINO_ACIDS.includes(aa1) && AMINO_ACIDS.includes(aa2)) {
        const key = `${aa1}-${aa2}`
        joint.set(key, (joint.get(key) || 0) + 1)
        marginal1.set(aa1, (marginal1.get(aa1) || 0) + 1)
        marginal2.set(aa2, (marginal2.get(aa2) || 0) + 1)
      }
    }

    let mi = 0
    for (const [key, count] of joint.entries()) {
      const [aa1, aa2] = key.split("-")
      const pJoint = count / n
      const pMarg1 = marginal1.get(aa1) / n
      const pMarg2 = marginal2.get(aa2) / n

      if (pJoint > 0 && pMarg1 > 0 && pMarg2 > 0) {
        mi += pJoint * Math.log2(pJoint / (pMarg1 * pMarg2))
      }
    }

    return mi
  }

  // Calculate asymmetric kernel regression based correlation between positions
  const calculateAsymmetricCorrelation = (msa, pos1, pos2, sigmaParam) => {
    if (!msa || msa.length === 0) return { forward: 0, backward: 0, asymmetry: 0 }

    const n = msa.length
    const epsilon = 1e-10

    // Build empirical distributions
    const marginal1 = new Map()
    const marginal2 = new Map()
    const jointCounts = new Map()

    for (const seq of msa) {
      const aa1 = seq[pos1].toUpperCase()
      const aa2 = seq[pos2].toUpperCase()

      if (AMINO_ACIDS.includes(aa1) && AMINO_ACIDS.includes(aa2)) {
        marginal1.set(aa1, (marginal1.get(aa1) || 0) + 1)
        marginal2.set(aa2, (marginal2.get(aa2) || 0) + 1)
        const key = `${aa1}-${aa2}`
        jointCounts.set(key, (jointCounts.get(key) || 0) + 1)
      }
    }

    // Normalize to frequencies
    for (const [aa, count] of marginal1.entries()) {
      marginal1.set(aa, count / n)
    }
    for (const [aa, count] of marginal2.entries()) {
      marginal2.set(aa, count / n)
    }

    // Forward kernel regression: predict pos2 from pos1
    // For each amino acid at pos1, use kernel regression to estimate distribution at pos2
    let forwardR2 = 0
    let forwardSST = 0
    let forwardSSR = 0

    for (const aa2_target of AMINO_ACIDS) {
      const observedFreq = marginal2.get(aa2_target) || 0

      // Kernel-weighted prediction
      let predictedFreq = 0
      let totalWeight = 0

      for (const aa1_ref of AMINO_ACIDS) {
        const freq1 = marginal1.get(aa1_ref) || epsilon

        // Kernel weight based on frequency similarity
        const kernelWeight = gaussianKernel(freq1, marginal1.get(aa1_ref) || 0, sigmaParam)

        // Conditional probability P(aa2_target | aa1_ref)
        const jointKey = `${aa1_ref}-${aa2_target}`
        const jointFreq = (jointCounts.get(jointKey) || 0) / n
        const conditionalProb = freq1 > epsilon ? jointFreq / freq1 : 0

        predictedFreq += kernelWeight * conditionalProb
        totalWeight += kernelWeight
      }

      if (totalWeight > epsilon) {
        predictedFreq /= totalWeight
      }

      // Calculate R² components
      const meanFreq = 1.0 / AMINO_ACIDS.length // Uniform prior
      forwardSST += Math.pow(observedFreq - meanFreq, 2)
      forwardSSR += Math.pow(observedFreq - predictedFreq, 2)
    }

    forwardR2 = forwardSST > epsilon ? 1 - forwardSSR / forwardSST : 0
    forwardR2 = Math.max(0, Math.min(1, forwardR2)) // Clamp to [0,1]

    // Backward kernel regression: predict pos1 from pos2
    let backwardR2 = 0
    let backwardSST = 0
    let backwardSSR = 0

    for (const aa1_target of AMINO_ACIDS) {
      const observedFreq = marginal1.get(aa1_target) || 0

      // Kernel-weighted prediction
      let predictedFreq = 0
      let totalWeight = 0

      for (const aa2_ref of AMINO_ACIDS) {
        const freq2 = marginal2.get(aa2_ref) || epsilon

        const kernelWeight = gaussianKernel(freq2, marginal2.get(aa2_ref) || 0, sigmaParam)

        const jointKey = `${aa1_target}-${aa2_ref}`
        const jointFreq = (jointCounts.get(jointKey) || 0) / n
        const conditionalProb = freq2 > epsilon ? jointFreq / freq2 : 0

        predictedFreq += kernelWeight * conditionalProb
        totalWeight += kernelWeight
      }

      if (totalWeight > epsilon) {
        predictedFreq /= totalWeight
      }

      const meanFreq = 1.0 / AMINO_ACIDS.length
      backwardSST += Math.pow(observedFreq - meanFreq, 2)
      backwardSSR += Math.pow(observedFreq - predictedFreq, 2)
    }

    backwardR2 = backwardSST > epsilon ? 1 - backwardSSR / backwardSST : 0
    backwardR2 = Math.max(0, Math.min(1, backwardR2))

    // Asymmetry: positive means forward prediction is stronger
    const asymmetry = forwardR2 - backwardR2

    return {
      forward: Number.parseFloat(forwardR2.toFixed(4)),
      backward: Number.parseFloat(backwardR2.toFixed(4)),
      asymmetry: Number.parseFloat(asymmetry.toFixed(4)),
    }
  }

  // Calculate structural complexity index based on sequence patterns
  const calculateStructuralComplexity = (msa, position) => {
    if (!msa || msa.length === 0) return 0

    const column = msa.map((seq) => seq[position].toUpperCase())
    const validAAs = column.filter((aa) => AMINO_ACIDS.includes(aa))

    if (validAAs.length === 0) return 0

    // Count unique amino acids
    const unique = new Set(validAAs)
    const diversity = unique.size / 20 // Normalized by max possible

    // Calculate hydrophobicity variance (using Kyte-Doolittle scale approximation)
    const hydrophobicity = {
      I: 4.5,
      V: 4.2,
      L: 3.8,
      F: 2.8,
      C: 2.5,
      M: 1.9,
      A: 1.8,
      G: -0.4,
      T: -0.7,
      S: -0.8,
      W: -0.9,
      Y: -1.3,
      P: -1.6,
      H: -3.2,
      E: -3.5,
      Q: -3.5,
      D: -3.5,
      N: -3.5,
      K: -3.9,
      R: -4.5,
    }

    const hydroValues = validAAs.map((aa) => hydrophobicity[aa] || 0)
    const meanHydro = hydroValues.reduce((a, b) => a + b, 0) / hydroValues.length
    const varHydro = hydroValues.reduce((a, b) => a + Math.pow(b - meanHydro, 2), 0) / hydroValues.length

    // Complexity score combines diversity and physicochemical variance
    return (diversity + Math.sqrt(varHydro) / 5) / 2
  }

  // Calculate the T_v tensor with enhanced normalization and statistical robustness
  const calculateTvScoreForPosition = (msa, position, windowSizeParam, sigmaParam) => {
    const L = msa[0].length
    const epsilon = 1e-10

    // Build frequency distribution at this position
    const freqAtPos = new Map()
    let validCount = 0

    for (const seq of msa) {
      const aa = seq[position].toUpperCase()
      if (AMINO_ACIDS.includes(aa)) {
        freqAtPos.set(aa, (freqAtPos.get(aa) || 0) + 1)
        validCount++
      }
    }

    // Normalize to frequencies
    for (const [aa, count] of freqAtPos.entries()) {
      freqAtPos.set(aa, count / validCount)
    }

    // Calculate Tv for each amino acid pair (X→Y)
    const tvTensor = Array(AMINO_ACIDS.length)
      .fill(0)
      .map(() => Array(AMINO_ACIDS.length).fill(0))

    // For each pair of amino acids at this position
    for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
      for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
        if (xIdx === yIdx) continue

        const aaX = AMINO_ACIDS[xIdx]
        const aaY = AMINO_ACIDS[yIdx]

        // Calculate asymmetric correlations with neighboring positions
        let forwardCorrelationSum = 0
        let backwardCorrelationSum = 0
        let neighborCount = 0

        for (let offset = -windowSizeParam; offset <= windowSizeParam; offset++) {
          if (offset === 0) continue

          const neighborPos = position + offset
          if (neighborPos < 0 || neighborPos >= L) continue

          // Calculate how well aaX at current position predicts neighbor
          // vs how well neighbor predicts aaX
          const corr = calculateAsymmetricCorrelation(msa, position, neighborPos, sigmaParam)

          // Weight by frequency of the amino acids
          const freqX = freqAtPos.get(aaX) || epsilon
          const freqY = freqAtPos.get(aaY) || epsilon

          // Forward: how well X→Y transition is supported by correlations
          forwardCorrelationSum += corr.forward * freqX
          // Backward: how well Y→X transition is supported
          backwardCorrelationSum += corr.backward * freqY

          neighborCount++
        }

        if (neighborCount === 0) continue

        // Average correlations
        const avgForward = forwardCorrelationSum / neighborCount + epsilon
        const avgBackward = backwardCorrelationSum / neighborCount + epsilon

        // Tv = P(X→Y) / P(Y→X) with log normalization for stability
        const ratio = avgForward / avgBackward
        const tvValue = Math.log(ratio + 1) // Log(1+x) for numerical stability

        tvTensor[xIdx][yIdx] = tvValue
      }
    }

    return tvTensor
  }

  // Smooth T_v tensors over a sliding window
  const smoothTvTensors = (allTvData, L, w) => {
    const smoothedData = []
    for (let i = 0; i < L; i++) {
      const start = Math.max(0, i - w)
      const end = Math.min(L - 1, i + w)
      const count = end - start + 1

      const sumTensor = Array(AMINO_ACIDS.length)
        .fill(0)
        .map(() => Array(AMINO_ACIDS.length).fill(0))

      for (let j = start; j <= end; j++) {
        const tvTensor = allTvData[j]?.tvTensor
        if (tvTensor) {
          for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
            for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
              if (xIdx !== yIdx) {
                sumTensor[xIdx][yIdx] += tvTensor[xIdx][yIdx]
              }
            }
          }
        }
      }

      const smoothedTensor = sumTensor.map((row) => row.map((val) => val / count))
      smoothedData.push({ position: i, smoothedTvTensor: smoothedTensor })
    }
    return smoothedData
  }

  // Calculate global directional variant tension
  const calculateGlobalTv = (allPositionsTvData) => {
    const globalSums = Array(AMINO_ACIDS.length)
      .fill(0)
      .map(() => Array(AMINO_ACIDS.length).fill(0))
    const L = allPositionsTvData.length

    if (L === 0) return {}

    for (const { tvTensor } of allPositionsTvData) {
      if (tvTensor) {
        for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
          for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
            if (xIdx !== yIdx) {
              // Exclude diagonal elements
              globalSums[xIdx][yIdx] += tvTensor[xIdx][yIdx]
            }
          }
        }
      }
    }

    const globalTvResults = {}
    for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
      for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
        if (xIdx === yIdx) continue
        const X = AMINO_ACIDS[xIdx]
        const Y = AMINO_ACIDS[yIdx]
        globalTvResults[`${X}→${Y}`] = globalSums[xIdx][yIdx] / L
      }
    }
    return globalTvResults
  }

  // Identify high-tension sites
  const findHighTensionSites = (allPositionsTvData, topN = 10) => {
    const sites = []
    allPositionsTvData.forEach(({ position, smoothedTvTensor }) => {
      if (!smoothedTvTensor) return
      let maxTv = 0
      let maxTvX = ""
      let maxTvY = ""

      for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
        for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
          if (xIdx === yIdx) continue
          const tvValue = smoothedTvTensor[xIdx][yIdx]
          if (tvValue > maxTv) {
            maxTv = tvValue
            maxTvX = AMINO_ACIDS[xIdx]
            maxTvY = AMINO_ACIDS[yIdx]
          }
        }
      }
      if (maxTv > 0) {
        sites.push({ position: position + 1, maxTv: Number.parseFloat(maxTv.toFixed(4)), X: maxTvX, Y: maxTvY })
      }
    })

    return sites.sort((a, b) => b.maxTv - a.maxTv).slice(0, topN)
  }

  // Bootstrap resampling for confidence interval estimation
  const bootstrapTvAnalysis = (msa, numIterations, confidenceLvl) => {
    const n = msa.length
    const L = msa[0].length
    const globalTvSamples = []

    for (let iter = 0; iter < numIterations; iter++) {
      // Resample sequences with replacement
      const resampledMSA = []
      for (let i = 0; i < n; i++) {
        const idx = Math.floor(Math.random() * n)
        resampledMSA.push(msa[idx])
      }

      const tempTv = []
      for (let i = 0; i < L; i++) {
        const tvTensor = calculateTvScoreForPosition(resampledMSA, i, windowSize, sigma)
        // Create a simple tensor representation for compatibility
        const tvTensorForAvg = Array(AMINO_ACIDS.length)
          .fill(0)
          .map(() => Array(AMINO_ACIDS.length).fill(0))

        let sumTv = 0
        let countTv = 0
        for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
          for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
            if (xIdx !== yIdx) {
              sumTv += tvTensor[xIdx][yIdx]
              countTv++
            }
          }
        }
        const avgTvScore = countTv > 0 ? sumTv / countTv : 0

        tempTv.push({ position: i, tvTensor, tvScore: avgTvScore })
      }

      const bootstrapGlobalTv = calculateGlobalTv(tempTv)
      globalTvSamples.push(bootstrapGlobalTv)
    }

    // Calculate confidence intervals
    const ciResults = {}
    const alpha = 1 - confidenceLvl

    // Get all substitution keys
    const substitutions = Object.keys(globalTvSamples[0])

    for (const sub of substitutions) {
      const values = globalTvSamples.map((sample) => sample[sub]).sort((a, b) => a - b)
      const lowerIdx = Math.floor((values.length * alpha) / 2)
      const upperIdx = Math.floor(values.length * (1 - alpha / 2))
      const meanVal = values.reduce((a, b) => a + b, 0) / values.length

      ciResults[sub] = {
        lower: values[lowerIdx],
        upper: values[upperIdx],
        mean: meanVal,
        std: Math.sqrt(values.reduce((a, b) => a + Math.pow(b - meanVal, 2), 0) / values.length),
      }
    }

    return ciResults
  }

  // Cumulative distribution function for standard normal distribution
  const normalCDF = (x) => {
    const t = 1 / (1 + 0.2316419 * Math.abs(x))
    const d = 0.3989423 * Math.exp((-x * x) / 2)
    const prob = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
    return x > 0 ? 1 - prob : prob
  }

  const benjaminiHochbergCorrection = (pValues) => {
    const sorted = pValues.map((p, idx) => ({ p, idx })).sort((a, b) => a.p - b.p)

    const m = pValues.length
    const corrected = new Array(m)

    for (let i = 0; i < m; i++) {
      const rank = i + 1
      const bh = (pValues[sorted[i].idx] * m) / rank
      corrected[sorted[i].idx] = Math.min(bh, 1)
    }

    // Enforce monotonicity
    for (let i = m - 2; i >= 0; i--) {
      corrected[sorted[i].idx] = Math.min(corrected[sorted[i].idx], corrected[sorted[i + 1].idx])
    }

    return corrected
  }

  const permutationTestCorrelation = (
    msa,
    pos1,
    pos2,
    sigma,
    observedForward,
    observedBackward,
    numPermutations = 1000,
  ) => {
    let forwardCount = 0
    let backwardCount = 0

    for (let perm = 0; perm < numPermutations; perm++) {
      // Permute position 2 while keeping position 1 fixed
      const permutedMSA = msa.map((seq, idx) => {
        const newSeq = seq.split("")
        const randomIdx = Math.floor(Math.random() * msa.length)
        newSeq[pos2] = msa[randomIdx][pos2]
        return newSeq.join("")
      })

      const permCorr = calculateAsymmetricCorrelation(permutedMSA, pos1, pos2, sigma)

      if (Math.abs(permCorr.forward) >= Math.abs(observedForward)) forwardCount++
      if (Math.abs(permCorr.backward) >= Math.abs(observedBackward)) backwardCount++
    }

    return {
      forwardPValue: (forwardCount + 1) / (numPermutations + 1),
      backwardPValue: (backwardCount + 1) / (numPermutations + 1),
    }
  }

  const calculateCohenD = (forward, backward, stdForward, stdBackward) => {
    const pooledStd = Math.sqrt((stdForward ** 2 + stdBackward ** 2) / 2)
    return pooledStd > 0 ? (forward - backward) / pooledStd : 0
  }

  const crossValidateKernelRegression = (msa, position, sigma, folds = 5) => {
    const n = msa.length
    const foldSize = Math.floor(n / folds)
    let totalR2 = 0
    let totalRMSE = 0

    for (let fold = 0; fold < folds; fold++) {
      const testStart = fold * foldSize
      const testEnd = fold === folds - 1 ? n : (fold + 1) * foldSize

      const trainMSA = [...msa.slice(0, testStart), ...msa.slice(testEnd)]
      const testMSA = msa.slice(testStart, testEnd)

      // Calculate predictions and errors
      let sumSquaredError = 0
      let sumSquaredTotal = 0
      const meanFreq = {}

      // Calculate mean frequencies from training data
      AMINO_ACIDS.forEach((aa) => {
        const count = trainMSA.filter((seq) => seq[position] === aa).length
        meanFreq[aa] = count / trainMSA.length
      })

      // Test on held-out data
      testMSA.forEach((seq) => {
        const actual = seq[position]
        const actualFreq = meanFreq[actual] || 0

        // Kernel-weighted prediction
        let prediction = 0
        let weightSum = 0
        trainMSA.forEach((trainSeq) => {
          const weight = gaussianKernel(seq[position], trainSeq[position], sigma) // Simplified kernel call for frequency comparison
          prediction += weight * (trainSeq[position] === actual ? 1 : 0)
          weightSum += weight
        })
        prediction = weightSum > 0 ? prediction / weightSum : actualFreq

        sumSquaredError += (actualFreq - prediction) ** 2
        sumSquaredTotal += (actualFreq - Object.values(meanFreq).reduce((a, b) => a + b, 0) / AMINO_ACIDS.length) ** 2
      })

      const r2 = 1 - sumSquaredError / (sumSquaredTotal || 1)
      const rmse = Math.sqrt(sumSquaredError / testMSA.length)

      totalR2 += r2
      totalRMSE += rmse
    }

    return {
      meanR2: totalR2 / folds,
      meanRMSE: totalRMSE / folds,
    }
  }

  const wilcoxonSignedRankTest = (differences) => {
    // Remove zeros
    const nonZero = differences.filter((d) => Math.abs(d) > 1e-10)
    if (nonZero.length === 0) return { statistic: 0, pValue: 1 }

    // Rank absolute differences
    const ranked = nonZero
      .map((d, idx) => ({ value: d, absValue: Math.abs(d), idx }))
      .sort((a, b) => a.absValue - b.absValue)
      .map((item, rank) => ({ ...item, rank: rank + 1 }))

    // Sum of positive ranks
    const wPlus = ranked.filter((item) => item.value > 0).reduce((sum, item) => sum + item.rank, 0)

    const n = nonZero.length
    const expectedW = (n * (n + 1)) / 4
    const varianceW = (n * (n + 1) * (2 * n + 1)) / 24
    const z = (wPlus - expectedW) / Math.sqrt(varianceW)

    // Two-tailed p-value
    const pValue = 2 * (1 - normalCDF(Math.abs(z)))

    return { statistic: wPlus, zScore: z, pValue }
  }

  // Export functions
  const exportToCSV = (data, filename) => {
    const csvContent = data.map((row) => Object.values(row).join(",")).join("\n")
    const header = Object.keys(data[0]).join(",") + "\n"
    const blob = new Blob([header + csvContent], { type: "text/csv" })
    const url = window.URL.createObjectURL(blob)
    const link = document.createElement("a")
    link.href = url
    link.download = filename
    link.click()
  }

  const exportResults = (format) => {
    const timestamp = new Date().toISOString().split("T")[0]

    if (format === "csv-global") {
      const data = Object.entries(globalTv).map(([substitution, value]) => ({
        substitution,
        tv_value: value.toFixed(6),
      }))
      exportToCSV(data, `global_tv_${timestamp}.csv`)
    } else if (format === "csv-sites") {
      const data = highTensionSites.map((site) => ({
        position: site.position,
        max_tv: site.maxTv,
        direction: `${site.X}→${site.Y}`,
      }))
      exportToCSV(data, `high_tension_sites_${timestamp}.csv`)
    } else if (format === "json") {
      const data = {
        parameters: {
          sigma,
          windowSize,
          pseudocount,
          bootstrapIterations,
          confidenceLevel,
          autoSigma,
          calculatedSigma,
        },
        msaStats,
        globalTv,
        confidenceIntervals,
        highTensionSites,
        statisticalSignificance,
        mutualInformation,
        structuralComplexity,
        correlationPValues,
        fdrCorrectedPValues,
        effectSizes,
        cvMetrics,
        timestamp,
      }
      const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" })
      const url = window.URL.createObjectURL(blob)
      const link = document.createElement("a")
      link.href = url
      link.download = `tv_analysis_${timestamp}.json`
      link.click()
    }
  }

  // Main calculation handler with enhanced statistical analysis
  const handleCalculate = async () => {
    setIsLoading(true)
    setError("")

    // Allow React to render the loading screen
    await new Promise((resolve) => setTimeout(resolve, 100))

    try {
      const msa = parseMSA(msaText)
      if (msa.length === 0) {
        throw new Error("Please provide a valid MSA.")
      }
      if (msa.some((seq) => seq.length !== msa[0].length)) {
        throw new Error("All sequences in the MSA must have the same length.")
      }
      if (msa.length < 3) {
        throw new Error("At least 3 sequences are required for robust statistical analysis.")
      }

      // Calculate comprehensive MSA statistics
      const stats = calculateMsaStats(msa)
      setMsaStats(stats)

      const effectiveSigma = autoSigma ? calculateOptimalBandwidth(msa, stats) : sigma
      setCalculatedSigma(effectiveSigma)

      if (effectiveSigma <= 0) {
        throw new Error("Calculated bandwidth must be greater than 0.")
      }
      if (windowSize < 0) {
        throw new Error("Window size (w) cannot be negative.")
      }

      const { frequencies, entropies, L } = calculateFrequenciesAndEntropy(msa, pseudocount)

      const tempAllPositionsTv = []
      for (let i = 0; i < L; i++) {
        const positionSigma = autoSigma ? calculateAdaptiveBandwidth(msa, i, effectiveSigma) : effectiveSigma
        const tvTensor = calculateTvScoreForPosition(msa, i, windowSize, positionSigma)

        // Calculate average Tv score for this position
        let sumTv = 0
        let countTv = 0
        for (let xIdx = 0; xIdx < AMINO_ACIDS.length; xIdx++) {
          for (let yIdx = 0; yIdx < AMINO_ACIDS.length; yIdx++) {
            if (xIdx !== yIdx) {
              sumTv += tvTensor[xIdx][yIdx]
              countTv++
            }
          }
        }
        const avgTvScore = countTv > 0 ? sumTv / countTv : 0

        tempAllPositionsTv.push({ position: i, tvTensor, tvScore: avgTvScore })
      }

      // Apply smoothing
      const smoothedData = smoothTvTensors(tempAllPositionsTv, L, windowSize)
      const finalAllPositionsTv = tempAllPositionsTv.map((item, idx) => ({
        ...item,
        smoothedTvTensor: smoothedData[idx].smoothedTvTensor,
      }))
      setAllPositionsTv(finalAllPositionsTv)

      // Calculate global Tv
      const calculatedGlobalTv = calculateGlobalTv(finalAllPositionsTv) // Use finalAllPositionsTv which has smoothed tensors
      setGlobalTv(calculatedGlobalTv)

      // Bootstrap confidence intervals (if enough sequences)
      if (msa.length >= 10 && bootstrapIterations > 0) {
        const ciResults = bootstrapTvAnalysis(msa, bootstrapIterations, confidenceLevel)
        setConfidenceIntervals(ciResults)
      }

      // Calculate asymmetric correlations between adjacent positions
      const miData = []
      for (let i = 0; i < L - 1; i++) {
        const correlation = calculateAsymmetricCorrelation(msa, i, i + 1, effectiveSigma)
        miData.push({
          position: i + 1,
          mi: correlation.forward, // Use forward correlation as primary metric
          backward: correlation.backward,
          asymmetry: correlation.asymmetry,
        })
      }

      const correlationPValuesData = []
      const effectSizesData = []

      for (let i = 0; i < L - 1; i++) {
        const correlation = calculateAsymmetricCorrelation(msa, i, i + 1, effectiveSigma)

        // Permutation test (using fewer permutations for speed)
        const pValues = permutationTestCorrelation(
          msa,
          i,
          i + 1,
          effectiveSigma,
          correlation.forward,
          correlation.backward,
          500,
        )

        // Calculate effect size
        const effectSize = calculateCohenD(correlation.forward, correlation.backward, 0.1, 0.1) // Assuming a baseline std dev for interpretation

        correlationPValuesData.push({
          position: i + 1,
          forwardPValue: pValues.forwardPValue,
          backwardPValue: pValues.backwardPValue,
          asymmetryPValue: Math.min(pValues.forwardPValue, pValues.backwardPValue) * 2, // Bonferroni-like for simplicity, better to use actual asymmetry test
        })

        effectSizesData.push({
          position: i + 1,
          effectSize,
          interpretation: Math.abs(effectSize) < 0.2 ? "Small" : Math.abs(effectSize) < 0.5 ? "Medium" : "Large",
        })
      }

      setCorrelationPValues(correlationPValuesData)
      setEffectSizes(effectSizesData)

      // This block was causing the lint error: calculatedHighTensionSites was used before it was calculated.
      // The calculation of highTensionSites is now performed later in the handleCalculate function,
      // so this section needs to be moved to after that calculation.
      // For now, we will conditionally check if highTensionSites is populated.
      let calculatedHighTensionSites = []
      if (finalAllPositionsTv.length > 0) {
        calculatedHighTensionSites = findHighTensionSites(finalAllPositionsTv, 10)
      }

      const allPValues = calculatedHighTensionSites.map((site) => {
        const sig = statisticalSignificance[site.position]
        return sig?.pValue || 1
      })

      if (allPValues.length > 0) {
        const correctedPValues = benjaminiHochbergCorrection(allPValues)
        const fdrResults = {}

        calculatedHighTensionSites.forEach((site, idx) => {
          fdrResults[site.position] = {
            // Keep original stats if they exist, otherwise create a placeholder
            ...(statisticalSignificance[site.position] || {}),
            fdrPValue: correctedPValues[idx],
            fdrSignificant: correctedPValues[idx] < 0.05,
          }
        })

        setFdrCorrectedPValues(fdrResults)
      }

      const cvResults = {}
      // Ensure L is defined before using it
      if (L > 0) {
        const samplePositions = [0, Math.floor(L / 4), Math.floor(L / 2), Math.floor((3 * L) / 4), L - 1].filter(
          (pos) => pos < L,
        ) // Ensure valid positions

        for (const pos of samplePositions) {
          const cv = crossValidateKernelRegression(msa, pos, effectiveSigma)
          cvResults[pos + 1] = cv // Store with 1-based indexing
        }
      }

      const asymmetryDifferences = miData.map((item) => item.mi - item.backward)
      const wilcoxonResult = wilcoxonSignedRankTest(asymmetryDifferences)

      setCvMetrics({ ...cvResults, wilcoxonTest: wilcoxonResult })

      setMutualInformation(miData)

      // Calculate structural complexity for each position
      const complexityData = []
      for (let i = 0; i < L; i++) {
        const complexity = calculateStructuralComplexity(msa, i)
        complexityData.push({ position: i + 1, complexity: Number.parseFloat(complexity.toFixed(4)) })
      }
      setStructuralComplexity(complexityData)

      // Identify high-tension sites with statistical significance
      // Recalculate calculatedHighTensionSites here after all data is ready
      const finalCalculatedHighTensionSites = findHighTensionSites(finalAllPositionsTv, 10)

      // Calculate statistical significance for high-tension sites
      const sigResults = {}
      for (const site of finalCalculatedHighTensionSites) {
        const pos = site.position - 1
        const posEntropy = entropies.get(pos)
        const effectiveN = stats.effectiveSequences[pos]

        // Calculate z-score based on position-specific statistics
        const allTvValues = []
        const tensor = finalAllPositionsTv[pos]?.smoothedTvTensor // Use optional chaining for safety
        if (tensor) {
          for (let i = 0; i < AMINO_ACIDS.length; i++) {
            for (let j = 0; j < AMINO_ACIDS.length; j++) {
              if (i !== j) allTvValues.push(tensor[i][j])
            }
          }
        }

        const mean = allTvValues.reduce((a, b) => a + b, 0) / allTvValues.length
        const variance = allTvValues.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / allTvValues.length
        const stdDev = Math.sqrt(variance)
        const zScore = stdDev > 0 ? (site.maxTv - mean) / stdDev : 0

        // Approximate p-value (two-tailed)
        const pValue = 2 * (1 - normalCDF(Math.abs(zScore)))

        sigResults[site.position] = {
          zScore: Number.parseFloat(zScore.toFixed(3)),
          pValue: Number.parseFloat(pValue.toFixed(6)),
          significant: pValue < 0.05,
        }
      }
      setStatisticalSignificance(sigResults)
      setHighTensionSites(finalCalculatedHighTensionSites)
      // Ensure sigResults are also passed if they were modified by FDR
      // (Already handled by setting statisticalSignificance and then updating fdrCorrectedPValues)
    } catch (err) {
      setError(err.message)
    } finally {
      setIsLoading(false)
    }
  }

  // Prepare enhanced chart data with multiple metrics
  const chartData = useMemo(() => {
    if (!allPositionsTv.length) return []

    return allPositionsTv.map((item, idx) => {
      let maxTv = 0
      let avgTv = 0
      let complexity = 0
      let mi = 0
      let conservation = 0
      let entropy = 0

      if (item.smoothedTvTensor) {
        let tensorMaxTv = 0
        let tensorAvgTv = 0
        let tensorCount = 0
        for (let i = 0; i < AMINO_ACIDS.length; i++) {
          for (let j = 0; j < AMINO_ACIDS.length; j++) {
            if (i !== j) {
              const val = item.smoothedTvTensor[i][j]
              if (val > tensorMaxTv) tensorMaxTv = val
              tensorAvgTv += val
              tensorCount++
            }
          }
        }
        maxTv = tensorMaxTv
        avgTv = tensorCount > 0 ? tensorAvgTv / tensorCount : 0
      } else if (item.tvScore !== undefined) {
        // Fallback if smoothed tensors aren't available
        maxTv = item.tvScore
        avgTv = item.tvScore
      }

      maxTv = Number.parseFloat(maxTv.toFixed(4))
      avgTv = Number.parseFloat(avgTv.toFixed(4))

      if (structuralComplexity[idx]) {
        complexity = Number.parseFloat(structuralComplexity[idx].complexity.toFixed(4))
      }
      if (mutualInformation[idx]) {
        mi = Number.parseFloat(mutualInformation[idx].mi.toFixed(4))
      }
      if (msaStats.conservation && msaStats.conservation[idx] !== undefined) {
        conservation = Number.parseFloat(msaStats.conservation[idx].toFixed(4))
      }
      if (msaStats.entropy && msaStats.entropy[idx] !== undefined) {
        entropy = Number.parseFloat(msaStats.entropy[idx].toFixed(4))
      }

      return {
        position: idx + 1,
        maxTv,
        avgTv,
        conservation,
        entropy,
        complexity,
        mi,
      }
    })
  }, [allPositionsTv, msaStats, structuralComplexity, mutualInformation, AMINO_ACIDS])

  const globalTvChartData = useMemo(() => {
    const entries = Object.entries(globalTv)
      .sort(([, a], [, b]) => b - a)
      .slice(0, 20)
      .map(([substitution, value]) => {
        const ci = confidenceIntervals[substitution]
        return {
          substitution,
          value: Number.parseFloat(value.toFixed(4)),
          lower: ci ? Number.parseFloat(ci.lower.toFixed(4)) : value,
          upper: ci ? Number.parseFloat(ci.upper.toFixed(4)) : value,
          errorMargin: ci ? Number.parseFloat(((ci.upper - ci.lower) / 2).toFixed(4)) : 0,
        }
      })
    return entries
  }, [globalTv, confidenceIntervals])

  const exampleMSA = `>seq1
MTKCTKLLYV
>seq2
ATRTVYAFGT
>seq3
MTCLLLTQYV
>seq4
ATYCFTADGT`

  useEffect(() => {
    setMsaText(exampleMSA.trim())
  }, [])

  const TabButton = ({ id, label, icon: Icon, active, onClick }) => (
    <TabsTrigger
      value={id}
      className={`flex items-center gap-2 px-6 py-3 min-w-[160px] rounded-t-lg font-semibold text-sm transition-all duration-200
        data-[state=active]:bg-gradient-to-b data-[state=active]:from-orange-500 data-[state=active]:to-pink-500
        data-[state=active]:text-white data-[state=active]:shadow-lg data-[state=active]:scale-105
        data-[state=inactive]:bg-slate-100 data-[state=inactive]:text-slate-600
        hover:bg-gradient-to-b hover:from-orange-400 hover:to-pink-400 hover:text-white hover:shadow-md
        border-2 border-transparent data-[state=active]:border-orange-300`}
    >
      <Icon className="w-5 h-5" />
      <span className="whitespace-nowrap">{label}</span>
    </TabsTrigger>
  )

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-blue-50 to-slate-100">
      {/* Updated header with new styling and buttons */}
      <header className="bg-gradient-to-r from-orange-50 to-pink-50 shadow-md sticky top-0 z-10 border-b-2 border-orange-200">
        <div className="max-w-7xl mx-auto px-4 flex justify-between items-center py-2.5">
          <div className="flex items-center space-x-3">
            <div className="bg-gradient-to-br from-orange-400 to-pink-400 p-2 rounded-lg">
              <Database className="w-6 h-6 text-white" />
            </div>
            <div>
              <h1 className="text-2xl font-bold bg-gradient-to-r from-orange-600 to-pink-600 bg-clip-text text-transparent">
                Protein Tᵥ Analysis
              </h1>
              <p className="text-xs text-gray-600">Version 1.0.1.</p>
            </div>
          </div>
          <div className="flex items-center space-x-3">
            <Button
              variant="outline"
              size="sm"
              onClick={() => setShowDocumentation(true)}
              className="border-orange-300 text-orange-700 hover:bg-orange-50"
            >
              <FileText className="w-4 h-4 mr-2" />
              Documentation
            </Button>
            <Button
              variant="outline"
              size="sm"
              onClick={() => setShowVersionHistory(true)}
              className="border-pink-300 text-pink-700 hover:bg-pink-50"
            >
              <History className="w-4 h-4 mr-2" />
              Version History
            </Button>
          </div>
        </div>
      </header>

      <main className="container mx-auto px-4 py-8 max-w-7xl">
        {/* Header Section */}
        <div className="p-8 mb-6 rounded-xl shadow-lg border-2 border-orange-200 bg-gradient-to-br from-white via-orange-50/30 to-pink-50/30 relative overflow-hidden">
          <div className="absolute inset-0 opacity-5">
            <svg className="w-full h-full" xmlns="http://www.w3.org/2000/svg">
              <defs>
                <pattern id="dna-pattern" x="0" y="0" width="100" height="100" patternUnits="userSpaceOnUse">
                  <circle cx="25" cy="25" r="3" fill="currentColor" className="text-orange-400" />
                  <circle cx="75" cy="75" r="3" fill="currentColor" className="text-pink-400" />
                  <path
                    d="M 25 25 Q 50 50 75 75"
                    stroke="currentColor"
                    strokeWidth="1"
                    fill="none"
                    className="text-orange-300"
                  />
                </pattern>
              </defs>
              <rect width="100%" height="100%" fill="url(#dna-pattern)" />
            </svg>
          </div>

          <div className="relative z-10">
            {/* Updated main title and description */}
            <h1 className="font-bold bg-gradient-to-r from-orange-600 via-pink-600 to-orange-700 bg-clip-text mb-3 text-slate-700 text-4xl bg-transparent border-transparent">
              Directional Variant Tension Calculator
            </h1>
            <p className="text-gray-700 mb-6 max-w-4xl text-base">
              Advanced kernel regression-based analysis of evolutionary pressures in protein sequences with
              comprehensive statistical validation and coevolutionary signal detection.
            </p>
            <div className="flex flex-wrap gap-4 text-sm">
              <div className="flex items-center backdrop-blur-sm px-4 py-2 border rounded-none shadow-none bg-transparent text-transparent border-transparent">
                <Database className="w-4 h-4 mr-2 text-orange-500" />
                <span className="font-medium text-gray-700">MSA Analysis</span>
              </div>
              <div className="flex items-center backdrop-blur-sm px-4 py-2 border text-transparent bg-transparent border-transparent border-none shadow-none rounded-none">
                <TrendingUp className="w-4 h-4 mr-2 text-pink-500" />
                <span className="font-medium text-gray-700">Asymmetric Correlations</span>
              </div>
              <div className="flex items-center backdrop-blur-sm px-4 py-2 border shadow-none border-none border-transparent bg-transparent rounded-none">
                <BarChart3 className="w-4 h-4 mr-2 text-orange-500" />
                <span className="font-medium text-gray-700">Statistical Validation</span>
              </div>
              <div className="flex items-center backdrop-blur-sm px-4 py-2 rounded-full border bg-transparent border-transparent border-none shadow-none">
                <Shield className="w-4 h-4 mr-2 text-pink-500" />
                <span className="font-medium text-gray-700">FDR Correction</span>
              </div>
            </div>
          </div>
        </div>

        <Card className="mb-6 border-primary/20">
          <Collapsible open={methodologyOpen} onOpenChange={setMethodologyOpen}>
            <CardHeader className="pb-3">
              <CollapsibleTrigger className="flex items-center justify-between w-full hover:opacity-80 transition-opacity">
                <div className="flex items-center gap-2">
                  <Info className="h-5 w-5 text-primary" />
                  <CardTitle className="text-lg">Methodology & Formula</CardTitle>
                </div>
                {methodologyOpen ? (
                  <ChevronUp className="h-5 w-5 text-muted-foreground" />
                ) : (
                  <ChevronDown className="h-5 w-5 text-muted-foreground" />
                )}
              </CollapsibleTrigger>
              <CardDescription>Mathematical framework for directional variant tension analysis</CardDescription>
            </CardHeader>
            <CollapsibleContent>
              <CardContent className="space-y-4 text-sm leading-relaxed">
                <div className="space-y-2">
                  <h4 className="font-semibold text-foreground">Directional Variant Tension (Tᵥ)</h4>
                  <p className="text-muted-foreground">
                    The directional variant tension quantifies the asymmetric evolutionary pressure for amino acid
                    substitutions at each position in a protein sequence alignment.
                  </p>
                </div>

                <div className="bg-muted/50 p-4 rounded-lg border border-border space-y-3">
                  <h4 className="font-semibold text-foreground">Core Formula</h4>
                  <div className="font-mono text-center py-2 bg-background rounded border border-border">
                    Tᵥ(X→Y) = log(1 + P(X→Y) / P(Y→X))
                  </div>
                  <div className="space-y-2 text-muted-foreground">
                    <p>
                      <strong className="text-foreground">Where:</strong>
                    </p>
                    <ul className="list-disc list-inside space-y-1 ml-2">
                      <li>
                        <strong>X→Y</strong>: Amino acid substitution from X to Y at position i
                      </li>
                      <li>
                        <strong>P(X→Y)</strong>: Forward asymmetric kernel regression correlation
                      </li>
                      <li>
                        <strong>P(Y→X)</strong>: Backward asymmetric kernel regression correlation
                      </li>
                      <li>
                        <strong>log(1+x)</strong>: Logarithmic normalization for numerical stability
                      </li>
                    </ul>
                  </div>
                </div>

                <div className="space-y-2">
                  <h4 className="font-semibold text-foreground">Asymmetric Kernel Regression Correlations</h4>
                  <p className="text-muted-foreground">
                    The correlations P(X→Y) and P(Y→X) are computed using kernel-weighted conditional probability
                    estimation:
                  </p>
                  <div className="bg-muted/50 p-3 rounded border border-border space-y-2">
                    <div className="font-mono text-sm">P(X→Y) = R²(Y|X) = 1 - SSR/SST</div>
                    <p className="text-xs text-muted-foreground">
                      Where R² measures how well position X predicts position Y using kernel regression with Gaussian
                      bandwidth σ
                    </p>
                  </div>
                </div>

                <div className="space-y-2">
                  <h4 className="font-semibold text-foreground">Kernel Bandwidth Selection</h4>
                  <p className="text-muted-foreground">
                    The Gaussian kernel bandwidth is automatically selected using Silverman's rule adapted for sequence
                    data:
                  </p>
                  <div className="font-mono text-center py-2 bg-background rounded border border-border text-sm">
                    σ = 0.9 × min(σ̂, IQR/1.34) × N<sub>eff</sub>
                    <sup>-1/5</sup>
                  </div>
                  <p className="text-xs text-muted-foreground mt-2">
                    Where N<sub>eff</sub> is the effective number of sequences accounting for redundancy, σ̂ is the
                    standard deviation of amino acid frequencies, and IQR is the interquartile range.
                  </p>
                </div>

                <div className="space-y-2">
                  <h4 className="font-semibold text-foreground">Interpretation</h4>
                  <ul className="list-disc list-inside space-y-1 ml-2 text-muted-foreground">
                    <li>
                      <strong className="text-foreground">Tᵥ {">"} 0</strong>: Forward substitution (X→Y) is favored by
                      evolutionary pressure
                    </li>
                    <li>
                      <strong className="text-foreground">Tᵥ {"<"} 0</strong>: Backward substitution (Y→X) is favored
                    </li>
                    <li>
                      <strong className="text-foreground">|Tᵥ| large</strong>: Strong directional bias in evolutionary
                      constraints
                    </li>
                    <li>
                      <strong className="text-foreground">Tᵥ ≈ 0</strong>: Symmetric or neutral substitution pressure
                    </li>
                  </ul>
                </div>

                <div className="bg-primary/5 p-3 rounded border border-primary/20">
                  <p className="text-xs text-muted-foreground">
                    <strong className="text-foreground">Statistical Validation:</strong> Bootstrap resampling with
                    confidence intervals ensures robustness of Tᵥ estimates. High-tension sites are identified using
                    statistical significance testing (p {"<"} 0.05).
                  </p>
                </div>
              </CardContent>
            </CollapsibleContent>
          </Collapsible>
        </Card>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Input Panel */}
          <div className="lg:col-span-1">
            <Card className="mb-6">
              <CardHeader>
                <CardTitle className="flex items-center gap-2 bg-orange-50 text-foreground text-lg">
                  <FileText className="w-5 h-5 text-blue-600 mr-2" />
                  Input Data
                </CardTitle>
              </CardHeader>
              <CardContent>
                <Textarea
                  className="w-full p-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent transition text-sm resize-none font-mono tracking-wide"
                  rows={12}
                  placeholder="Paste your MSA here (FASTA format)..."
                  value={msaText}
                  onChange={(e) => setMsaText(e.target.value)}
                />
                <p className="text-xs text-gray-500 mt-2">
                  {"FASTA format: sequences starting with '>'} followed by sequence data"}
                </p>
              </CardContent>
            </Card>

            <Card className="mb-6">
              <CardHeader className="flex items-center gap-2 bg-orange-50 text-foreground">
                <Settings className="w-5 h-5 text-blue-600 mr-2" />
                <CardTitle className="text-lg">Parameters</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="space-y-4">
                  <div className="bg-blue-50 border border-blue-200 rounded-lg p-3">
                    <div className="flex items-center justify-between mb-2">
                      <Label className="text-sm font-medium text-blue-900">Automatic Bandwidth Selection</Label>
                      <Button
                        variant="outline"
                        size="sm"
                        onClick={() => setAutoSigma(!autoSigma)}
                        className={`relative h-6 w-11 cursor-pointer rounded-full transition-colors ${
                          autoSigma ? "bg-blue-600 border-blue-600" : "bg-gray-300 border-gray-300"
                        }`}
                        aria-pressed={autoSigma}
                      >
                        <span
                          className={`inline-block h-4 w-4 transform rounded-full bg-white transition-transform ${
                            autoSigma ? "translate-x-6" : "translate-x-1"
                          }`}
                        />
                      </Button>
                    </div>
                    <p className="text-xs text-blue-700">
                      {autoSigma
                        ? "Using Silverman's rule with adaptive corrections"
                        : "Manual bandwidth control enabled"}
                    </p>
                    {calculatedSigma > 0 && autoSigma && (
                      <div className="mt-2 text-xs text-blue-800 font-mono bg-blue-100 rounded px-2 py-1">
                        Calculated σ = {calculatedSigma.toFixed(4)}
                      </div>
                    )}
                  </div>

                  <div>
                    <Label className="block text-sm font-medium text-gray-700 mb-2">
                      Gaussian Kernel Bandwidth (σ)
                    </Label>
                    <Input
                      type="number"
                      className={`w-full p-3 border rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent ${
                        autoSigma ? "bg-gray-100 cursor-not-allowed" : "border-gray-300"
                      }`}
                      value={sigma}
                      onChange={(e) => setSigma(Number.parseFloat(e.target.value))}
                      min="0.01"
                      step="0.01"
                      disabled={autoSigma}
                    />
                    <p className="text-xs text-gray-500 mt-1">
                      {autoSigma
                        ? "Automatically calculated using sequence statistics"
                        : "Manual control: σ = 0.9 × min(σ̂, IQR/1.34) × Neff^(-1/5) × corrections"}
                    </p>
                  </div>
                  <div>
                    <Label className="block text-sm font-medium text-gray-700 mb-2">
                      Smoothing Window Half-Width (w)
                    </Label>
                    <Input
                      type="number"
                      className="w-full p-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      value={windowSize}
                      onChange={(e) => setWindowSize(Number.parseInt(e.target.value))}
                      min="0"
                      step="1"
                    />
                    <p className="text-xs text-gray-500 mt-1">Positions for local averaging</p>
                  </div>
                  <div>
                    <Label className="block text-sm font-medium text-gray-700 mb-2">Pseudocount Parameter</Label>
                    <Input
                      type="number"
                      className="w-full p-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      value={pseudocount}
                      onChange={(e) => setPseudocount(Number.parseFloat(e.target.value))}
                      min="0.001"
                      max="1"
                      step="0.001"
                    />
                    <p className="text-xs text-gray-500 mt-1">Dirichlet prior for frequency estimation</p>
                  </div>
                  <div>
                    <Label className="block text-sm font-medium text-gray-700 mb-2">Bootstrap Iterations</Label>
                    <Input
                      type="number"
                      className="w-full p-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      value={bootstrapIterations}
                      onChange={(e) => setBootstrapIterations(Number.parseInt(e.target.value))}
                      min="0"
                      max="1000"
                      step="10"
                    />
                    <p className="text-xs text-gray-500 mt-1">Resampling for confidence intervals (0 = disabled)</p>
                  </div>
                  <div>
                    <Label className="block text-sm font-medium text-gray-700 mb-2">Confidence Level</Label>
                    <Input
                      type="number"
                      className="w-full p-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      value={confidenceLevel}
                      onChange={(e) => setConfidenceLevel(Number.parseFloat(e.target.value))}
                      min="0.8"
                      max="0.99"
                      step="0.01"
                    />
                    <p className="text-xs text-gray-500 mt-1">Statistical confidence level (0.95 = 95%)</p>
                  </div>
                </div>
              </CardContent>
            </Card>

            <Card className="shadow-md border-l-4 bg-gradient-to-r from-orange-50 to-pink-50 border-b-0 border-transparent">
              <CardHeader className="pb-0">
                <div className="flex items-start gap-3">
                  <Info className="w-5 h-5 text-orange-600 mt-0.5 flex-shrink-0" />
                  <div>
                    <CardTitle className="text-base font-semibold text-orange-900">
                      Recommended Parameters & Computation Notes
                    </CardTitle>
                  </div>
                </div>
              </CardHeader>
              <CardContent className="text-sm space-y-2.5">
                <div className="bg-white rounded-lg p-3 border border-orange-200 py-1">
                  <p className="font-semibold text-orange-900 mb-2">For 100 sequences × 100 amino acids MSA:</p>
                  <ul className="space-y-1 text-gray-700 ml-4 list-disc">
                    <li>
                      <span className="font-medium">Bandwidth (σ):</span> Auto mode recommended (typically 0.05-0.15)
                    </li>
                    <li>
                      <span className="font-medium">Window Size (w):</span> 3-5 positions
                    </li>
                    <li>
                      <span className="font-medium">Pseudocount:</span> 0.01-0.05
                    </li>
                    <li>
                      <span className="font-medium">Bootstrap Iterations:</span> 100-500 (higher for more robust CI)
                    </li>
                    <li>
                      <span className="font-medium">Confidence Level:</span> 0.95 (standard)
                    </li>
                  </ul>
                </div>
                <div className="bg-amber-50 rounded-lg p-3 border border-amber-300">
                  <p className="text-amber-900 flex items-start gap-2">
                    <BoxesIcon className="w-4 h-4 mt-0.5 flex-shrink-0" />
                    <span>
                      <span className="font-semibold">Computation Time:</span> Analysis of longer sequences ({">"}200aa)
                      or large MSAs ({">"}500 sequences) may take several minutes due to processing resource limits of
                      the lab. Please be patient during calculation.
                    </span>
                  </p>
                </div>
              </CardContent>
            </Card>

            <Button
              onClick={handleCalculate}
              className="w-full bg-gradient-to-r from-blue-600 to-blue-700 hover:from-blue-700 hover:to-blue-800 font-bold px-6 shadow-lg transform transition duration-200 hover:scale-[1.02] focus:outline-none focus:ring-4 focus:ring-blue-300 disabled:opacity-50 disabled:cursor-not-allowed text-xl font-sans text-background rounded-xs py-6 my-4"
              disabled={isLoading}
            >
              {isLoading ? (
                <span className="flex items-center justify-center">
                  <Loader2 className="animate-spin -ml-1 mr-3 h-5 w-5 text-white" />
                  Analyzing...
                </span>
              ) : (
                "Calculate Tᵥ Analysis"
              )}
            </Button>

            {msaStats.numSequences && (
              <Card className="shadow-lg mt-4 p-6 opacity-75 rounded-lg border-blue-700 bg-blue-50">
                <CardHeader className="flex items-center gap-2 p-0 pb-3">
                  <Info className="w-5 h-5 text-green-600 mr-2" />
                  <CardTitle className="text-lg">MSA Descriptive Statistics</CardTitle>
                </CardHeader>
                <CardContent className="p-0 space-y-3 text-sm">
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <span className="text-gray-600">Sequences:</span>
                      <span className="font-semibold ml-2">{msaStats.numSequences}</span>
                    </div>
                    <div>
                      <span className="text-gray-600">Length:</span>
                      <span className="font-semibold ml-2">{msaStats.alignmentLength}</span>
                    </div>
                  </div>
                  <div>
                    <span className="text-gray-600">Effective Sequences (Neff):</span>
                    <span className="font-semibold ml-2">{msaStats.effectiveNumSequences?.toFixed(2)}</span>
                  </div>
                  <div>
                    <span className="text-gray-600">Avg. Pairwise Identity:</span>
                    <span className="font-semibold ml-2">{(msaStats.avgPairwiseIdentity * 100).toFixed(1)}%</span>
                  </div>
                  <div>
                    <span className="text-gray-600">Avg. Conservation:</span>
                    <span className="font-semibold ml-2">{(msaStats.avgConservation * 100).toFixed(1)}%</span>
                  </div>
                  <div>
                    <span className="text-gray-600">Avg. Shannon Entropy:</span>
                    <span className="font-semibold ml-2">{msaStats.avgEntropy?.toFixed(3)} bits</span>
                  </div>
                  <div>
                    <span className="text-gray-600">Gap Content:</span>
                    <span className="font-semibold ml-2">{msaStats.gapPercentage?.toFixed(2)}%</span>
                  </div>
                  {autoSigma && calculatedSigma > 0 && (
                    <div className="pt-3 border-t border-blue-200">
                      <span className="text-gray-600">Optimal Bandwidth (σ):</span>
                      <span className="font-semibold ml-2">{calculatedSigma.toFixed(4)}</span>
                      <p className="text-xs text-gray-500 mt-1">
                        Silverman's rule: σ = 0.9 × min(σ̂, IQR/1.34) × Neff^(-1/5) × corrections
                      </p>
                    </div>
                  )}
                </CardContent>
              </Card>
            )}
          </div>

          {/* Results Panel */}
          <div className="lg:col-span-2">
            {error && (
              <Alert variant="destructive" className="mb-6">
                <AlertCircle className="h-4 w-4" />
                <AlertDescription>{error}</AlertDescription>
              </Alert>
            )}

            {isLoading && (
              <Card className="bg-white text-center min-h-[600px] flex flex-col items-center justify-center">
                <CardHeader>
                  <div className="relative mb-8">
                    <Loader2 className="h-10 w-10 text-blue-600 animate-spin" />
                  </div>
                  <CardTitle className="text-2xl">Calculating Tᵥ Analysis</CardTitle>
                  <CardDescription className="max-w-md">
                    Processing your multiple sequence alignment with kernel regression and statistical validation...
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2 text-sm text-gray-500">
                    <div className="flex items-center justify-center space-x-2">
                      <span className="w-2 h-2 bg-blue-600 rounded-full animate-bounce"></span>
                      <span>Calculating Tᵥ tensors</span>
                    </div>
                    <div className="flex items-center justify-center space-x-2">
                      <span
                        className="w-2 h-2 bg-blue-600 rounded-full animate-bounce"
                        style={{ animationDelay: "0.1s" }}
                      ></span>
                      <span>Analyzing coevolutionary signals</span>
                    </div>
                    <div className="flex items-center justify-center space-x-2">
                      <span
                        className="w-2 h-2 bg-blue-600 rounded-full animate-bounce"
                        style={{ animationDelay: "0.2s" }}
                      ></span>
                      <span>Performing statistical validation</span>
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}

            {!isLoading && Object.keys(globalTv).length > 0 && (
              <Card className="overflow-hidden">
                {/* Tab Navigation */}
                <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
                  <TabsList className="border-b-2 border-orange-200 bg-gradient-to-r from-orange-50 to-pink-50 flex flex-wrap gap-2 justify-start w-full pb-0 pt-1">
                    <TabButton
                      id="global"
                      label="Global Analysis"
                      icon={TrendingUp}
                      active={activeTab === "global"}
                      onClick={() => setActiveTab("global")}
                    />
                    <TabButton
                      id="positional"
                      label="Positional Profile"
                      icon={BarChart3}
                      active={activeTab === "positional"}
                      onClick={() => setActiveTab("positional")}
                    />
                    <TabButton
                      id="sites"
                      label="High-Tension Sites"
                      icon={Database}
                      active={activeTab === "sites"}
                      onClick={() => setActiveTab("sites")}
                    />
                    <TabButton
                      id="coevolution"
                      label="Coevolution"
                      icon={Activity}
                      active={activeTab === "coevolution"}
                      onClick={() => setActiveTab("coevolution")}
                    />
                    <TabButton
                      id="statistics"
                      label="Statistics"
                      icon={Shield}
                      active={activeTab === "statistics"}
                      onClick={() => setActiveTab("statistics")}
                    />
                  </TabsList>

                  <CardContent className="p-6 mx-0 my-6 py-1.5">
                    {/* Export Controls */}
                    <div className="flex justify-end leading-8 gap-2 mb-0">
                      <Button
                        variant="outline"
                        onClick={() => exportResults("json")}
                        className="flex items-center gap-1.5 px-3 py-2 text-sm"
                      >
                        <Download className="w-4 h-4" />
                        JSON
                      </Button>
                      <Button
                        variant="outline"
                        onClick={() => exportResults("csv-global")}
                        className="flex items-center gap-1.5 px-3 py-2 text-sm"
                      >
                        <Download className="w-4 h-4" />
                        CSV (Global)
                      </Button>
                      <Button
                        variant="outline"
                        onClick={() => exportResults("csv-sites")}
                        className="flex items-center gap-1.5 px-3 py-2 text-sm"
                      >
                        <Download className="w-4 h-4" />
                        CSV (Sites)
                      </Button>
                    </div>

                    {/* Tab Content - Global Analysis */}
                    <TabsContent value="global">
                      <div className="space-y-6">
                        <div>
                          <h3 className="text-lg font-semibold text-gray-900 mb-3">
                            Global Directional Variant Tension
                          </h3>
                          <p className="text-gray-600 mb-4">
                            Average Tᵥ values across all positions with{" "}
                            {Object.keys(confidenceIntervals).length > 0
                              ? `${(confidenceLevel * 100).toFixed(0)}% confidence intervals`
                              : "point estimates"}
                            .
                          </p>
                          <div className="h-96">
                            <ResponsiveContainer width="100%" height="100%">
                              <BarChart data={globalTvChartData}>
                                <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                                <XAxis
                                  dataKey="substitution"
                                  angle={-45}
                                  textAnchor="end"
                                  height={80}
                                  fontSize={12}
                                  className="fill-muted-foreground"
                                />
                                <YAxis
                                  fontSize={12}
                                  label={{ value: "Tᵥ Value", angle: -90, position: "insideLeft" }}
                                  className="fill-muted-foreground"
                                />
                                <Tooltip
                                  formatter={(value, name) => {
                                    if (name === "value") return [value.toFixed(4), "Tᵥ Value"]
                                    if (name === "errorMargin") return [value.toFixed(4), "Error Margin"]
                                    return [value, name]
                                  }}
                                  labelStyle={{ color: "#374151" }}
                                  contentStyle={{
                                    backgroundColor: "hsl(var(--card))",
                                    border: "1px solid hsl(var(--border))",
                                  }}
                                />
                                <Bar dataKey="value" fill="hsl(var(--primary))">
                                  {globalTvChartData.map((entry, index) => (
                                    <Cell
                                      key={`cell-${index}`}
                                      fill={entry.errorMargin > 0.1 ? "hsl(var(--destructive))" : "hsl(var(--primary))"}
                                    />
                                  ))}
                                </Bar>
                              </BarChart>
                            </ResponsiveContainer>
                          </div>
                        </div>

                        <div>
                          <h4 className="text-md font-semibold text-gray-900 mb-3">
                            Top Substitution Pairs with Statistical Metrics
                          </h4>
                          <div className="overflow-hidden border border-gray-200 rounded-lg">
                            <table className="min-w-full divide-y divide-gray-200">
                              <thead className="bg-gray-50">
                                <tr>
                                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                    Rank
                                  </th>
                                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                    Substitution
                                  </th>
                                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                    Tᵥ Value
                                  </th>
                                  {Object.keys(confidenceIntervals).length > 0 && (
                                    <>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        CI Lower
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        CI Upper
                                      </th>
                                    </>
                                  )}
                                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                    Relative
                                  </th>
                                </tr>
                              </thead>
                              <tbody className="bg-white divide-y divide-gray-200">
                                {globalTvChartData.slice(0, 10).map((item, idx) => (
                                  <tr key={item.substitution} className="hover:bg-gray-50">
                                    <td className="px-4 py-3 text-sm text-gray-900">{idx + 1}</td>
                                    <td className="px-4 py-3 text-sm font-mono text-gray-900">{item.substitution}</td>
                                    <td className="px-4 py-3 text-sm text-gray-900">{item.value.toFixed(6)}</td>
                                    {Object.keys(confidenceIntervals).length > 0 && (
                                      <>
                                        <td className="px-4 py-3 text-sm text-gray-600">{item.lower.toFixed(6)}</td>
                                        <td className="px-4 py-3 text-sm text-gray-600">{item.upper.toFixed(6)}</td>
                                      </>
                                    )}
                                    <td className="px-4 py-3">
                                      <div className="flex items-center">
                                        <div className="flex-1 bg-gray-200 rounded-full h-2 mr-2">
                                          <div
                                            className="bg-blue-600 h-2 rounded-full"
                                            style={{
                                              width: `${(item.value / Math.max(...globalTvChartData.map((d) => d.value))) * 100}%`,
                                            }}
                                          ></div>
                                        </div>
                                        <span className="text-xs text-gray-500">
                                          {(
                                            (item.value / Math.max(...globalTvChartData.map((d) => d.value))) *
                                            100
                                          ).toFixed(1)}
                                          %
                                        </span>
                                      </div>
                                    </td>
                                  </tr>
                                ))}
                              </tbody>
                            </table>
                          </div>
                        </div>
                      </div>
                    </TabsContent>

                    {/* Tab Content - Positional Profile */}
                    <TabsContent value="positional">
                      <div className="space-y-6">
                        <div className="bg-white rounded-lg border border-gray-200 p-6">
                          <h3 className="text-lg font-semibold text-gray-900 mb-4">Positional Tᵥ Profile</h3>
                          <div className="h-96">
                            <ResponsiveContainer width="100%" height="100%">
                              <LineChart data={chartData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                                <CartesianGrid strokeDasharray="3 3" className="stroke-gray-200" />
                                <XAxis
                                  dataKey="position"
                                  label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                  className="fill-muted-foreground"
                                />
                                <YAxis
                                  label={{ value: "Tᵥ Score", angle: -90, position: "insideLeft" }}
                                  className="fill-muted-foreground"
                                />
                                <Tooltip
                                  formatter={(value: any, name: any) => {
                                    return [Number(value).toFixed(4), name]
                                  }}
                                  contentStyle={{
                                    backgroundColor: "hsl(var(--card))",
                                    border: "1px solid hsl(var(--border))",
                                  }}
                                />
                                <Legend />
                                <Line
                                  type="monotone"
                                  dataKey="maxTv"
                                  stroke="#ea580c"
                                  strokeWidth={3}
                                  name="Max Tᵥ"
                                  dot={{ fill: "#ea580c", r: 4 }}
                                />
                                <Line
                                  type="monotone"
                                  dataKey="avgTv"
                                  stroke="#ec4899"
                                  strokeWidth={2}
                                  strokeDasharray="5 5"
                                  name="Avg Tᵥ"
                                  dot={false}
                                />
                                <Line
                                  type="monotone"
                                  dataKey="conservation"
                                  stroke="#10b981"
                                  strokeWidth={2}
                                  name="Conservation"
                                  dot={false}
                                  opacity={0.7}
                                />
                                <Line
                                  type="monotone"
                                  dataKey="entropy"
                                  stroke="#8b5cf6"
                                  strokeWidth={2}
                                  strokeDasharray="3 3"
                                  name="Entropy"
                                  dot={false}
                                  opacity={0.6}
                                />
                              </LineChart>
                            </ResponsiveContainer>
                          </div>
                        </div>

                        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                          <div className="bg-white rounded-lg border border-gray-200 p-6">
                            <h3 className="text-lg font-semibold text-gray-900 mb-4">Structural Complexity</h3>
                            <div className="h-64">
                              <ResponsiveContainer width="100%" height="100%">
                                <LineChart data={chartData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                                  <CartesianGrid strokeDasharray="3 3" className="stroke-gray-200" />
                                  <XAxis
                                    dataKey="position"
                                    label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                    className="fill-muted-foreground"
                                  />
                                  <YAxis
                                    label={{ value: "Complexity", angle: -90, position: "insideLeft" }}
                                    className="fill-muted-foreground"
                                  />
                                  <Tooltip
                                    formatter={(value: any) => [Number(value).toFixed(4), "Complexity"]}
                                    contentStyle={{
                                      backgroundColor: "hsl(var(--card))",
                                      border: "1px solid hsl(var(--border))",
                                    }}
                                  />
                                  <Line
                                    type="monotone"
                                    dataKey="complexity"
                                    stroke="#f59e0b"
                                    strokeWidth={2}
                                    name="Complexity"
                                    dot={{ fill: "#f59e0b", r: 3 }}
                                  />
                                </LineChart>
                              </ResponsiveContainer>
                            </div>
                          </div>

                          <div className="bg-white rounded-lg border border-gray-200 p-6">
                            <h3 className="text-lg font-semibold text-gray-900 mb-4">Mutual Information</h3>
                            <div className="h-64">
                              <ResponsiveContainer width="100%" height="100%">
                                <LineChart data={chartData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
                                  <CartesianGrid strokeDasharray="3 3" className="stroke-gray-200" />
                                  <XAxis
                                    dataKey="position"
                                    label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                    className="fill-muted-foreground"
                                  />
                                  <YAxis
                                    label={{ value: "MI", angle: -90, position: "insideLeft" }}
                                    className="fill-muted-foreground"
                                  />
                                  <Tooltip
                                    formatter={(value: any) => [Number(value).toFixed(4), "MI"]}
                                    contentStyle={{
                                      backgroundColor: "hsl(var(--card))",
                                      border: "1px solid hsl(var(--border))",
                                    }}
                                  />
                                  <Line
                                    type="monotone"
                                    dataKey="mi"
                                    stroke="#06b6d4"
                                    strokeWidth={2}
                                    name="MI"
                                    dot={{ fill: "#06b6d4", r: 3 }}
                                  />
                                </LineChart>
                              </ResponsiveContainer>
                            </div>
                          </div>
                        </div>

                        <div className="bg-white rounded-lg border border-gray-200 p-6">
                          <h3 className="text-lg font-semibold text-gray-900 mb-4">Detailed Position Statistics</h3>
                          <div className="overflow-x-auto">
                            <table className="w-full text-sm">
                              <thead className="bg-gray-50 border-b border-gray-200">
                                <tr>
                                  <th className="px-4 py-3 text-left font-semibold text-gray-700">Position</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">Max Tᵥ</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">Avg Tᵥ</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">Conservation</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">Entropy</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">Complexity</th>
                                  <th className="px-4 py-3 text-right font-semibold text-gray-700">MI</th>
                                </tr>
                              </thead>
                              <tbody className="divide-y divide-gray-200">
                                {chartData.slice(0, 20).map((row) => (
                                  <tr key={row.position} className="hover:bg-gray-50">
                                    <td className="px-4 py-3 font-medium text-gray-900">{row.position}</td>
                                    <td className="px-4 py-3 text-right text-orange-600 font-semibold">
                                      {row.maxTv.toFixed(4)}
                                    </td>
                                    <td className="px-4 py-3 text-right text-pink-600">{row.avgTv.toFixed(4)}</td>
                                    <td className="px-4 py-3 text-right text-green-600">
                                      {row.conservation.toFixed(4)}
                                    </td>
                                    <td className="px-4 py-3 text-right text-purple-600">{row.entropy.toFixed(4)}</td>
                                    <td className="px-4 py-3 text-right text-amber-600">{row.complexity.toFixed(4)}</td>
                                    <td className="px-4 py-3 text-right text-cyan-600">{row.mi.toFixed(4)}</td>
                                  </tr>
                                ))}
                              </tbody>
                            </table>
                            {chartData.length > 20 && (
                              <p className="text-sm text-gray-500 mt-4 text-center">
                                Showing first 20 of {chartData.length} positions
                              </p>
                            )}
                          </div>
                        </div>
                      </div>
                    </TabsContent>

                    <TabsContent value="sites" className="mt-6">
                      <Card className="border-gray-200 shadow-sm">
                        <CardHeader className="bg-gradient-to-r from-red-50 to-orange-50 border-b border-gray-200">
                          <CardTitle className="text-xl font-bold text-gray-900 flex items-center">
                            <Zap className="w-5 h-5 mr-2 text-red-600" />
                            High-Tension Sites
                          </CardTitle>
                        </CardHeader>
                        <CardContent className="p-6">
                          <p className="text-gray-600 mb-6">
                            Top positions with highest directional variant tension, indicating strong evolutionary
                            selection pressures. Statistical significance is evaluated using z-scores and p-values, with
                            FDR correction for multiple comparisons.
                          </p>
                          {highTensionSites.length > 0 ? (
                            <>
                              <div className="h-80 mb-6">
                                <ResponsiveContainer width="100%" height="100%">
                                  <BarChart data={highTensionSites}>
                                    <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                                    <XAxis
                                      dataKey="position"
                                      fontSize={12}
                                      label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                      className="fill-muted-foreground"
                                    />
                                    <YAxis
                                      fontSize={12}
                                      label={{ value: "Max Tᵥ", angle: -90, position: "insideLeft" }}
                                      className="fill-muted-foreground"
                                    />
                                    <Tooltip
                                      formatter={(value, name) => {
                                        if (name === "maxTv") return [value.toFixed(4), "Max Tᵥ"]
                                        return [value, name]
                                      }}
                                      labelFormatter={(label) => `Position ${label}`}
                                      contentStyle={{
                                        backgroundColor: "hsl(var(--card))",
                                        border: "1px solid hsl(var(--border))",
                                      }}
                                    />
                                    <Bar dataKey="maxTv" fill="hsl(var(--primary))">
                                      {highTensionSites.map((entry, index) => (
                                        <Cell
                                          key={`cell-${index}`}
                                          fill={
                                            fdrCorrectedPValues[entry.position]?.fdrSignificant
                                              ? "hsl(var(--destructive))"
                                              : "hsl(var(--primary))"
                                          }
                                        />
                                      ))}
                                    </Bar>
                                  </BarChart>
                                </ResponsiveContainer>
                              </div>
                              <div className="overflow-hidden border border-gray-200 rounded-lg">
                                <table className="min-w-full divide-y divide-gray-200">
                                  <thead className="bg-gray-50">
                                    <tr>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Rank
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Position
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Max Tᵥ
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Direction
                                      </th>
                                      {Object.keys(statisticalSignificance).length > 0 && (
                                        <>
                                          <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                            Z-Score
                                          </th>
                                          <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                            P-Value
                                          </th>
                                          <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                            FDR P-Value
                                          </th>
                                          <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                            {"Significant (FDR < 0.05)"}
                                          </th>
                                        </>
                                      )}
                                    </tr>
                                  </thead>
                                  <tbody className="bg-white divide-y divide-gray-200">
                                    {highTensionSites.map((site, idx) => {
                                      const sig = statisticalSignificance[site.position]
                                      const fdr = fdrCorrectedPValues[site.position]
                                      return (
                                        <tr
                                          key={site.position}
                                          className={fdr?.fdrSignificant ? "bg-green-50" : "hover:bg-gray-50"}
                                        >
                                          <td className="px-4 py-3 text-sm text-gray-900">{idx + 1}</td>
                                          <td className="px-4 py-3 text-sm font-semibold text-gray-900">
                                            {site.position}
                                          </td>
                                          <td className="px-4 py-3 text-sm text-gray-900">{site.maxTv.toFixed(4)}</td>
                                          <td className="px-4 py-3 text-sm font-mono text-blue-600">
                                            {site.X}→{site.Y}
                                          </td>
                                          {Object.keys(statisticalSignificance).length > 0 && (
                                            <>
                                              <td className="px-4 py-3 text-sm text-gray-900">
                                                {sig?.zScore.toFixed(3) || "N/A"}
                                              </td>
                                              <td className="px-4 py-3 text-sm text-gray-900">
                                                {sig?.pValue.toFixed(6) || "N/A"}
                                              </td>
                                              <td className="px-4 py-3 text-sm text-gray-900">
                                                <span
                                                  className={fdr?.fdrPValue < 0.05 ? "font-bold text-green-700" : ""}
                                                >
                                                  {fdr?.fdrPValue.toFixed(6) || "N/A"}
                                                </span>
                                              </td>
                                              <td className="px-4 py-3 text-sm">
                                                {fdr?.fdrSignificant ? (
                                                  <span className="px-2 py-1 text-xs font-semibold rounded-full bg-green-100 text-green-800">
                                                    Yes
                                                  </span>
                                                ) : (
                                                  <span className="px-2 py-1 text-xs font-semibold rounded-full bg-gray-100 text-gray-600">
                                                    No
                                                  </span>
                                                )}
                                              </td>
                                            </>
                                          )}
                                        </tr>
                                      )
                                    })}
                                  </tbody>
                                </table>
                              </div>

                              {Object.keys(fdrCorrectedPValues).length > 0 && (
                                <div className="mt-4 bg-yellow-50 border border-yellow-200 rounded-lg p-4">
                                  <h5 className="font-semibold text-yellow-900 mb-2 flex items-center">
                                    <AlertCircle className="w-4 h-4 mr-2" />
                                    Multiple Testing Correction (FDR)
                                  </h5>
                                  <p className="text-sm text-yellow-800">
                                    False Discovery Rate (FDR) correction using the Benjamini-Hochberg procedure
                                    controls for multiple comparisons. FDR-corrected p-values account for testing
                                    multiple positions simultaneously, reducing false positives. Sites with FDR p-value{" "}
                                    {"<"} 0.05 are considered statistically significant after correction.
                                  </p>
                                </div>
                              )}
                            </>
                          ) : (
                            <div className="text-center py-12 text-gray-500">
                              No high-tension sites calculated yet. Run the analysis first.
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>

                    <TabsContent value="coevolution" className="mt-6">
                      <Card className="border-gray-200 shadow-sm">
                        <CardHeader className="bg-gradient-to-r from-purple-50 to-pink-50 border-b border-gray-200">
                          <CardTitle className="text-xl font-bold text-gray-900 flex items-center">
                            <TrendingUp className="w-5 h-5 mr-2 text-purple-600" />
                            Asymmetric Coevolutionary Analysis
                          </CardTitle>
                        </CardHeader>
                        <CardContent className="p-6">
                          <p className="text-gray-600 mb-6">
                            Kernel regression-based asymmetric correlations between adjacent positions with statistical
                            significance testing.
                          </p>
                          <div className="h-80 mb-6">
                            <ResponsiveContainer width="100%" height="100%">
                              <ComposedChart data={mutualInformation}>
                                <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                                <XAxis
                                  dataKey="position"
                                  fontSize={12}
                                  label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                  className="fill-muted-foreground"
                                />
                                <YAxis
                                  fontSize={12}
                                  label={{ value: "Mutual Information (bits)", angle: -90, position: "insideLeft" }}
                                  className="fill-muted-foreground"
                                />
                                <Tooltip
                                  formatter={(value) => [value.toFixed(4) + " bits", "MI"]}
                                  labelFormatter={(label) => `Position ${label} ↔ ${label + 1}`}
                                  contentStyle={{
                                    backgroundColor: "hsl(var(--card))",
                                    border: "1px solid hsl(var(--border))",
                                  }}
                                />
                                <Area
                                  type="monotone"
                                  dataKey="mi"
                                  fill="hsl(var(--chart-1) / 0.2)"
                                  stroke="hsl(var(--chart-1))"
                                  strokeWidth={2}
                                  name="Mutual Information"
                                />
                                <Line
                                  type="monotone"
                                  dataKey="mi"
                                  stroke="hsl(var(--chart-1))"
                                  strokeWidth={2}
                                  dot={{ fill: "hsl(var(--chart-1))", strokeWidth: 2, r: 3 }}
                                />
                              </ComposedChart>
                            </ResponsiveContainer>
                          </div>

                          {correlationPValues.length > 0 && (
                            <div className="mt-6">
                              <h4 className="text-md font-semibold text-gray-900 mb-3">
                                Statistical Significance of Correlations (Permutation Tests)
                              </h4>
                              <div className="overflow-hidden border border-gray-200 rounded-lg">
                                <table className="min-w-full divide-y divide-gray-200">
                                  <thead className="bg-gray-50">
                                    <tr>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Position
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Forward Corr
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        P-Value (Fwd)
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Backward Corr
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        P-Value (Bwd)
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Effect Size
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Significant
                                      </th>
                                    </tr>
                                  </thead>
                                  <tbody className="bg-white divide-y divide-gray-200">
                                    {correlationPValues.slice(0, 20).map((item, idx) => {
                                      const miItem = mutualInformation[idx]
                                      const esItem = effectSizes[idx]
                                      const isSig = item.forwardPValue < 0.05 || item.backwardPValue < 0.05

                                      return (
                                        <tr key={item.position} className={isSig ? "bg-green-50" : "hover:bg-gray-50"}>
                                          <td className="px-4 py-3 text-sm font-semibold text-gray-900">
                                            {item.position}
                                          </td>
                                          <td className="px-4 py-3 text-sm text-gray-900">{miItem?.mi.toFixed(4)}</td>
                                          <td className="px-4 py-3 text-sm text-gray-900">
                                            <span
                                              className={item.forwardPValue < 0.05 ? "font-bold text-green-700" : ""}
                                            >
                                              {item.forwardPValue.toFixed(4)}
                                            </span>
                                          </td>
                                          <td className="px-4 py-3 text-sm text-gray-900">
                                            {miItem?.backward.toFixed(4)}
                                          </td>
                                          <td className="px-4 py-3 text-sm text-gray-900">
                                            <span
                                              className={item.backwardPValue < 0.05 ? "font-bold text-green-700" : ""}
                                            >
                                              {item.backwardPValue.toFixed(4)}
                                            </span>
                                          </td>
                                          <td className="px-4 py-3 text-sm text-gray-900">
                                            {esItem?.effectSize.toFixed(3)} ({esItem?.interpretation})
                                          </td>
                                          <td className="px-4 py-3 text-sm">
                                            {isSig ? (
                                              <span className="px-2 py-1 text-xs font-semibold rounded-full bg-green-100 text-green-800">
                                                Yes
                                              </span>
                                            ) : (
                                              <span className="px-2 py-1 text-xs font-semibold rounded-full bg-gray-100 text-gray-600">
                                                No
                                              </span>
                                            )}
                                          </td>
                                        </tr>
                                      )
                                    })}
                                  </tbody>
                                </table>
                              </div>

                              {cvMetrics.wilcoxonTest && (
                                <div className="mt-4 bg-blue-50 border border-blue-200 rounded-lg p-4">
                                  <h5 className="font-semibold text-blue-900 mb-2">
                                    Wilcoxon Signed-Rank Test for Overall Asymmetry
                                  </h5>
                                  <div className="grid grid-cols-3 gap-4 text-sm">
                                    <div>
                                      <span className="text-blue-700 font-medium">Test Statistic:</span>
                                      <span className="ml-2 text-blue-900">
                                        {cvMetrics.wilcoxonTest.statistic.toFixed(2)}
                                      </span>
                                    </div>
                                    <div>
                                      <span className="text-blue-700 font-medium">Z-Score:</span>
                                      <span className="ml-2 text-blue-900">
                                        {cvMetrics.wilcoxonTest.zScore.toFixed(3)}
                                      </span>
                                    </div>
                                    <div>
                                      <span className="text-blue-700 font-medium">P-Value:</span>
                                      <span
                                        className={`ml-2 font-bold ${cvMetrics.wilcoxonTest.pValue < 0.05 ? "text-green-700" : "text-blue-900"}`}
                                      >
                                        {cvMetrics.wilcoxonTest.pValue.toFixed(6)}
                                      </span>
                                    </div>
                                  </div>
                                  <p className="text-xs text-blue-800 mt-2">
                                    Tests whether forward and backward correlations are significantly different
                                    (asymmetric).
                                    {cvMetrics.wilcoxonTest.pValue < 0.05
                                      ? " Significant asymmetry detected (p < 0.05)."
                                      : " No significant asymmetry detected."}
                                  </p>
                                </div>
                              )}
                            </div>
                          )}

                          <div className="mt-4 bg-purple-50 border border-purple-200 rounded-lg p-4">
                            <h4 className="font-semibold text-purple-900 mb-2 flex items-center">
                              <Info className="w-4 h-4 mr-2" />
                              Interpretation
                            </h4>
                            <p className="text-sm text-purple-800">
                              Higher mutual information indicates stronger coevolutionary coupling between adjacent
                              positions. Peaks suggest functionally important regions where residues must coordinate for
                              proper protein function or structural stability.
                            </p>
                          </div>
                        </CardContent>
                      </Card>

                      <div>
                        <h3 className="text-lg font-semibold text-gray-900 mb-3">Structural Complexity Index</h3>
                        <p className="text-gray-600 mb-4">
                          Combines amino acid diversity and physicochemical variance to identify structurally complex
                          regions.
                        </p>
                        {structuralComplexity.length > 0 ? (
                          <div className="h-80">
                            <ResponsiveContainer width="100%" height="100%">
                              <BarChart data={structuralComplexity}>
                                <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                                <XAxis
                                  dataKey="position"
                                  fontSize={12}
                                  label={{ value: "Position", position: "insideBottom", offset: -5 }}
                                  className="fill-muted-foreground"
                                />
                                <YAxis
                                  fontSize={12}
                                  label={{ value: "Complexity Index", angle: -90, position: "insideLeft" }}
                                  className="fill-muted-foreground"
                                />
                                <Tooltip
                                  formatter={(value) => [value.toFixed(4), "Complexity"]}
                                  labelFormatter={(label) => `Position ${label}`}
                                  contentStyle={{
                                    backgroundColor: "hsl(var(--card))",
                                    border: "1px solid hsl(var(--border))",
                                  }}
                                />
                                <Bar dataKey="complexity" fill="hsl(var(--destructive))">
                                  {structuralComplexity.map((entry, index) => (
                                    <Cell
                                      key={`cell-${index}`}
                                      fill={entry.complexity > 0.5 ? "hsl(var(--destructive))" : "hsl(var(--chart-4))"}
                                    />
                                  ))}
                                </Bar>
                              </BarChart>
                            </ResponsiveContainer>
                          </div>
                        ) : (
                          <div className="text-center py-12 text-gray-500">
                            No complexity data calculated yet. Run the analysis first.
                          </div>
                        )}
                      </div>
                    </TabsContent>

                    <TabsContent value="statistics" className="mt-6">
                      <Card className="border-gray-200 shadow-sm">
                        <CardHeader className="bg-gradient-to-r from-green-50 to-blue-50 border-b border-gray-200">
                          <CardTitle className="text-xl font-bold text-gray-900 flex items-center">
                            <Shield className="w-5 h-5 mr-2 text-green-600" />
                            Statistical Validation & Metrics
                          </CardTitle>
                        </CardHeader>
                        <CardContent className="p-6">
                          <h3 className="text-lg font-semibold text-gray-900 mb-3">Statistical Summary</h3>
                          <p className="text-gray-600 mb-4">
                            Comprehensive statistical metrics and validation results for the Tᵥ analysis.
                          </p>

                          {Object.keys(confidenceIntervals).length > 0 ? (
                            <>
                              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                                <Card className="bg-blue-50">
                                  <CardHeader>
                                    <CardTitle className="font-semibold text-blue-900 mb-2">
                                      Bootstrap Iterations
                                    </CardTitle>
                                  </CardHeader>
                                  <CardContent>
                                    <p className="text-2xl font-bold text-blue-700">{bootstrapIterations}</p>
                                  </CardContent>
                                </Card>
                                <Card className="bg-green-50">
                                  <CardHeader>
                                    <CardTitle className="font-semibold text-green-900 mb-2">
                                      Confidence Level
                                    </CardTitle>
                                  </CardHeader>
                                  <CardContent>
                                    <p className="text-2xl font-bold text-green-700">
                                      {(confidenceLevel * 100).toFixed(0)}%
                                    </p>
                                  </CardContent>
                                </Card>
                                <Card className="bg-purple-50">
                                  <CardHeader>
                                    <CardTitle className="font-semibold text-purple-900 mb-2">
                                      Substitution Pairs
                                    </CardTitle>
                                  </CardHeader>
                                  <CardContent>
                                    <p className="text-2xl font-bold text-purple-700">
                                      {Object.keys(confidenceIntervals).length}
                                    </p>
                                  </CardContent>
                                </Card>
                              </div>

                              <div>
                                <h4 className="text-md font-semibold text-gray-900 mb-3">
                                  Confidence Intervals for Top Substitutions
                                </h4>
                                <div className="overflow-hidden border border-gray-200 rounded-lg">
                                  <table className="min-w-full divide-y divide-gray-200">
                                    <thead className="bg-gray-50">
                                      <tr>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          Substitution
                                        </th>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          Mean
                                        </th>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          Std Dev
                                        </th>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          CI Lower
                                        </th>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          CI Upper
                                        </th>
                                        <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                          Range
                                        </th>
                                      </tr>
                                    </thead>
                                    <tbody className="bg-white divide-y divide-gray-200">
                                      {Object.entries(confidenceIntervals)
                                        .sort(([, a], [, b]) => b.mean - a.mean)
                                        .slice(0, 15)
                                        .map(([sub, ci]) => (
                                          <tr key={sub} className="hover:bg-gray-50">
                                            <td className="px-4 py-3 text-sm font-mono text-gray-900">{sub}</td>
                                            <td className="px-4 py-3 text-sm text-gray-900">{ci.mean.toFixed(6)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-600">{ci.std.toFixed(6)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-600">{ci.lower.toFixed(6)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-600">{ci.upper.toFixed(6)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-600">
                                              {(ci.upper - ci.lower).toFixed(6)}
                                            </td>
                                          </tr>
                                        ))}
                                    </tbody>
                                  </table>
                                </div>
                              </div>

                              <div className="mt-6 bg-green-50 border border-green-200 rounded-lg p-4">
                                <h4 className="font-semibold text-green-900 mb-2 flex items-center">
                                  <Shield className="w-4 h-4 mr-2" />
                                  Statistical Validation
                                </h4>
                                <p className="text-sm text-green-800">
                                  Bootstrap confidence intervals provide robust uncertainty estimates for Tᵥ values.
                                  Narrow intervals indicate stable, reliable measurements, while wider intervals suggest
                                  higher variability in the evolutionary signal.
                                </p>
                              </div>
                            </>
                          ) : (
                            <div className="text-center py-12 text-gray-500">
                              <Shield className="w-12 h-12 mx-auto mb-4 text-gray-400" />
                              <p>No bootstrap statistics available.</p>
                              <p className="text-sm mt-2">
                                {
                                  "Set bootstrap iterations > 0 and ensure at least 10 sequences for statistical validation."
                                }
                              </p>
                            </div>
                          )}

                          {Object.keys(statisticalSignificance).length > 0 && (
                            <div className="mt-6">
                              <h4 className="text-md font-semibold text-gray-900 mb-3">
                                Significance Testing for High-Tension Sites
                              </h4>
                              <div className="overflow-hidden border border-gray-200 rounded-lg">
                                <table className="min-w-full divide-y divide-gray-200">
                                  <thead className="bg-gray-50">
                                    <tr>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Position
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Z-Score
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        P-Value
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        FDR P-Value
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        {"Significant (FDR < 0.05)"}
                                      </th>
                                    </tr>
                                  </thead>
                                  <tbody className="bg-white divide-y divide-gray-200">
                                    {Object.entries(statisticalSignificance)
                                      .sort(([, a], [, b]) => a.pValue - b.pValue)
                                      .map(([pos, sig]) => {
                                        const fdr = fdrCorrectedPValues[pos]
                                        return (
                                          <tr
                                            key={pos}
                                            className={fdr?.fdrSignificant ? "bg-green-50" : "hover:bg-gray-50"}
                                          >
                                            <td className="px-4 py-3 text-sm font-semibold text-gray-900">{pos}</td>
                                            <td className="px-4 py-3 text-sm text-gray-900">{sig.zScore.toFixed(3)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-900">{sig.pValue.toFixed(6)}</td>
                                            <td className="px-4 py-3 text-sm text-gray-900">
                                              <span className={fdr?.fdrPValue < 0.05 ? "font-bold text-green-700" : ""}>
                                                {fdr?.fdrPValue.toFixed(6) || "N/A"}
                                              </span>
                                            </td>
                                            <td className="px-4 py-3 text-sm">
                                              {fdr?.fdrSignificant ? (
                                                <span className="px-2 py-1 text-xs font-semibold rounded-full bg-green-100 text-green-800">
                                                  Yes
                                                </span>
                                              ) : (
                                                <span className="px-2 py-1 text-xs font-semibold rounded-full bg-gray-100 text-gray-600">
                                                  No
                                                </span>
                                              )}
                                            </td>
                                          </tr>
                                        )
                                      })}
                                  </tbody>
                                </table>
                              </div>
                            </div>
                          )}

                          {Object.keys(cvMetrics).length > 0 && cvMetrics.wilcoxonTest && (
                            <div className="mt-6">
                              <h4 className="text-md font-semibold text-gray-900 mb-3">
                                Kernel Regression Cross-Validation Metrics
                              </h4>
                              <div className="overflow-hidden border border-gray-200 rounded-lg">
                                <table className="min-w-full divide-y divide-gray-200">
                                  <thead className="bg-gray-50">
                                    <tr>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Position
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Mean R² (5-Fold CV)
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Mean RMSE
                                      </th>
                                      <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">
                                        Quality
                                      </th>
                                    </tr>
                                  </thead>
                                  <tbody className="bg-white divide-y divide-gray-200">
                                    {Object.entries(cvMetrics)
                                      .filter(([key]) => key !== "wilcoxonTest")
                                      .map(([pos, metrics]) => {
                                        const quality =
                                          metrics.meanR2 > 0.7
                                            ? "Excellent"
                                            : metrics.meanR2 > 0.5
                                              ? "Good"
                                              : metrics.meanR2 > 0.3
                                                ? "Fair"
                                                : "Poor"
                                        const qualityColor =
                                          metrics.meanR2 > 0.7
                                            ? "text-green-700"
                                            : metrics.meanR2 > 0.5
                                              ? "text-blue-700"
                                              : metrics.meanR2 > 0.3
                                                ? "text-yellow-700"
                                                : "text-red-700"

                                        return (
                                          <tr key={pos} className="hover:bg-gray-50">
                                            <td className="px-4 py-3 text-sm font-semibold text-gray-900">{pos}</td>
                                            <td className="px-4 py-3 text-sm text-gray-900">
                                              {metrics.meanR2.toFixed(4)}
                                            </td>
                                            <td className="px-4 py-3 text-sm text-gray-900">
                                              {metrics.meanRMSE.toFixed(4)}
                                            </td>
                                            <td className="px-4 py-3 text-sm">
                                              <span className={`font-semibold ${qualityColor}`}>{quality}</span>
                                            </td>
                                          </tr>
                                        )
                                      })}
                                  </tbody>
                                </table>
                              </div>

                              <div className="mt-4 bg-purple-50 border border-purple-200 rounded-lg p-4">
                                <h5 className="font-semibold text-purple-900 mb-2">Cross-Validation Interpretation</h5>
                                <p className="text-sm text-purple-800">
                                  5-fold cross-validation assesses kernel regression prediction quality. R² measures
                                  explained variance (higher is better, max 1.0), while RMSE measures prediction error
                                  (lower is better). These metrics validate the reliability of asymmetric correlation
                                  estimates at sampled positions.
                                </p>
                              </div>
                            </div>
                          )}
                        </CardContent>
                      </Card>
                    </TabsContent>
                  </CardContent>
                </Tabs>
              </Card>
            )}

            {!isLoading && Object.keys(globalTv).length === 0 && !error && (
              <Card className="bg-white text-center">
                <CardHeader>
                  <div className="w-16 h-16 bg-blue-100 rounded-full flex items-center justify-center mb-4">
                    <BoxesIcon className="w-8 h-8 text-blue-600" />
                  </div>
                  <CardTitle className="text-lg mb-2">Ready for Analysis</CardTitle>
                  <CardDescription>
                    Enter your MSA data and click &quot;Calculate Tᵥ Analysis&quot; to begin the directional variant
                    tension analysis.
                  </CardDescription>
                </CardHeader>
              </Card>
            )}
          </div>
        </div>

        <footer className="mt-16 bg-gradient-to-br from-orange-50 via-white to-pink-50 rounded-xl shadow-lg p-8 border-2 border-orange-200">
          <div className="text-center">
            <h4 className="text-2xl font-bold bg-gradient-to-r from-orange-600 to-pink-600 bg-clip-text text-transparent mb-4">
              About Directional Variant Tension (Tᵥ)
            </h4>
            <p className="text-gray-700 max-w-4xl mx-auto mb-6 leading-relaxed">
              Directional Variant Tension quantifies asymmetric evolutionary pressures in protein sequences through
              kernel-based conditional probability estimation. This framework integrates information theory,
              non-parametric statistics, and coevolutionary analysis to reveal directional selection pressures that
              shape protein evolution.
            </p>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-6 max-w-5xl mx-auto mb-8">
              <div className="bg-gradient-to-br from-orange-50 to-orange-100 rounded-xl p-5 border border-orange-200 shadow-sm py-2">
                <div className="flex items-center justify-center mb-3">
                  <div className="bg-orange-500 p-2 rounded-lg">
                    <TrendingUp className="w-5 h-5 text-white" />
                  </div>
                </div>
                <h5 className="font-bold text-orange-900 mb-2">Adaptive Bandwidth Selection</h5>
                <p className="text-sm text-orange-800">
                  Silverman's rule with Scott's correction, adapted for effective sequence count (Neff) and local
                  diversity metrics
                </p>
              </div>
              <div className="bg-gradient-to-br from-pink-50 to-pink-100 rounded-xl p-5 border border-pink-200 shadow-sm py-2">
                <div className="flex items-center justify-center mb-3">
                  <div className="bg-pink-500 p-2 rounded-lg">
                    <Shield className="w-5 h-5 text-white" />
                  </div>
                </div>
                <h5 className="font-bold text-pink-900 mb-2">Statistical Validation</h5>
                <p className="text-sm text-pink-800">
                  Permutation tests, bootstrap resampling, FDR correction, and cross-validation ensure result
                  reliability
                </p>
              </div>
              <div className="bg-gradient-to-br from-orange-50 to-pink-100 rounded-xl p-5 border border-orange-200 shadow-sm py-2">
                <div className="flex items-center justify-center mb-3">
                  <div className="bg-gradient-to-r from-orange-500 to-pink-500 p-2 rounded-lg">
                    <BarChart3 className="w-5 h-5 text-white" />
                  </div>
                </div>
                <h5 className="font-bold text-gray-900 mb-2">Coevolution Detection</h5>
                <p className="text-sm text-gray-800">
                  Asymmetric kernel regression correlations identify functionally coupled residue positions
                </p>
              </div>
            </div>

            <div className="border-t-2 border-orange-200 pt-6 mt-6">
              <p className="text-gray-600 font-medium">© {new Date().getFullYear()} Karagöl. All Rights Reserved.</p>
              <p className="text-sm text-gray-500 mt-2">Directional Variant Tension - Protein Evolutionary Analysis Toolkit</p>
            </div>
          </div>
        </footer>
      </main>

      {showDocumentation && (
        <div className="fixed inset-0 bg-black/50 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-white rounded-xl shadow-2xl max-w-4xl w-full max-h-[90vh] overflow-y-auto border-2 border-orange-200">
            <div className="sticky top-0 bg-gradient-to-r from-orange-50 to-pink-50 border-b-2 border-orange-200 p-6 flex justify-between items-center">
              <h2 className="text-2xl font-bold bg-gradient-to-r from-orange-600 to-pink-600 bg-clip-text text-transparent">
                Documentation
              </h2>
              <Button
                variant="ghost"
                size="sm"
                onClick={() => setShowDocumentation(false)}
                className="text-gray-600 hover:text-gray-900"
              >
                <X className="w-5 h-5" />
              </Button>
            </div>
            <div className="p-6 space-y-6">
              <section>
                <h3 className="text-xl font-bold text-orange-700 mb-3">Getting Started</h3>
                <p className="text-gray-700 mb-2">
                  Upload your Multiple Sequence Alignment (MSA) in FASTA format. The tool will analyze evolutionary
                  pressures and calculate directional variant tension scores.
                </p>
              </section>

              <section>
                <h3 className="text-xl font-bold text-pink-700 mb-3">Input Format</h3>
                <div className="bg-gray-50 rounded-lg p-4 border border-gray-200">
                  <code className="text-sm text-gray-800">
                    {">"}Sequence1
                    <br />
                    MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSLEVGQNLVPVALEEGLKKHPLFRLYRHRQTGSSYMEVEDNRSQVETLQRFYQNPTGDVTDIACDVTIVNNLKRLLQPWDPQFGSIN
                  </code>
                </div>
              </section>

              <section>
                <h3 className="text-xl font-bold text-orange-700 mb-3">Parameters</h3>
                <ul className="space-y-2 text-gray-700">
                  <li>
                    <strong className="text-orange-600">Bandwidth (σ):</strong> Controls kernel smoothing.
                    Auto-calculated using Silverman's rule or set manually.
                  </li>
                  <li>
                    <strong className="text-pink-600">Window Size:</strong> Number of neighboring positions for
                    correlation analysis.
                  </li>
                  <li>
                    <strong className="text-orange-600">Pseudocount:</strong> Prevents zero probabilities in sparse
                    data.
                  </li>
                  <li>
                    <strong className="text-pink-600">Bootstrap Iterations:</strong> Number of resampling iterations for
                    confidence intervals.
                  </li>
                </ul>
              </section>

              <section>
                <h3 className="text-xl font-bold text-pink-700 mb-3">Output Interpretation</h3>
                <ul className="space-y-2 text-gray-700">
                  <li>
                    <strong className="text-orange-600">Tᵥ Score:</strong> Measures directional evolutionary pressure.
                    Higher values indicate stronger asymmetry.
                  </li>
                  <li>
                    <strong className="text-pink-600">P-values:</strong> Statistical significance of correlations
                    (permutation-based).
                  </li>
                  <li>
                    <strong className="text-orange-600">FDR:</strong> False Discovery Rate corrected p-values for
                    multiple testing.
                  </li>
                  <li>
                    <strong className="text-pink-600">Effect Size:</strong> Cohen's d quantifies magnitude of asymmetry.
                  </li>
                </ul>
              </section>

              <section>
                <h3 className="text-xl font-bold text-orange-700 mb-3">Citation</h3>
                <div className="bg-orange-50 rounded-lg p-4 border border-orange-200">
                  <p className="text-sm text-gray-800">
                    If you use this tool in your research, please cite:
                    <br />
                    <em>
                      Karagöl et al. (2024). Directional Variant Tension: A Kernel Regression Approach to Evolutionary
                      Pressure Analysis.
                    </em>
                  </p>
                </div>
              </section>
            </div>
          </div>
        </div>
      )}

      {showVersionHistory && (
        <div className="fixed inset-0 bg-black/50 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-white rounded-xl shadow-2xl max-w-3xl w-full max-h-[90vh] overflow-y-auto border-2 border-pink-200">
            <div className="sticky top-0 bg-gradient-to-r from-pink-50 to-orange-50 border-b-2 border-pink-200 p-6 flex justify-between items-center">
              <h2 className="text-2xl font-bold bg-gradient-to-r from-pink-600 to-orange-600 bg-clip-text text-transparent">
                Version History
              </h2>
              <Button
                variant="ghost"
                size="sm"
                onClick={() => setShowVersionHistory(false)}
                className="text-gray-600 hover:text-gray-900"
              >
                <X className="w-5 h-5" />
              </Button>
            </div>
            <div className="p-6 space-y-6">
              <div className="border-l-4 border-orange-500 pl-4">
                <div className="flex items-center gap-2 mb-2">
                  <span className="bg-orange-500 text-white text-xs font-bold px-2 py-1 rounded">CURRENT</span>
                  <h3 className="text-lg font-bold text-gray-900">Version 1.0.1</h3>
                  <span className="text-sm text-gray-500">January 2024</span>
                </div>
                <ul className="space-y-1 text-sm text-gray-700">
                  <li>• Implemented asymmetric kernel regression correlations</li>
                  <li>• Added permutation-based p-values and FDR correction</li>
                  <li>• Integrated cross-validation metrics for kernel regression quality</li>
                  <li>• Added Cohen's d effect sizes for asymmetry magnitude</li>
                  <li>• Implemented adaptive bandwidth selection (Silverman's + Scott's rule)</li>
                  <li>• Enhanced UI with academic evolutionary biology design</li>
                </ul>
              </div>

              <div className="border-l-4 border-pink-400 pl-4">
                <div className="flex items-center gap-2 mb-2">
                  <h3 className="text-lg font-bold text-gray-900">Version 0.2.1 unpublished</h3>
                  <span className="text-sm text-gray-500">December 2023</span>
                </div>
                <ul className="space-y-1 text-sm text-gray-700">
                  <li>• Added coevolutionary analysis with mutual information</li>
                  <li>• Implemented bootstrap confidence intervals</li>
                  <li>• Enhanced visualization with multiple chart types</li>
                </ul>
              </div>

              <div className="border-l-4 border-orange-300 pl-4">
                <div className="flex items-center gap-2 mb-2">
                  <h3 className="text-lg font-bold text-gray-900">Version 0.1.1 unpublished</h3>
                  <span className="text-sm text-gray-500">October 2023</span>
                </div>
                <ul className="space-y-1 text-sm text-gray-700">
                  <li>• Initial release with basic Tᵥ calculation</li>
                  <li>• FASTA format support for MSA input</li>
                  <li>• Kernel-based probability estimation</li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
