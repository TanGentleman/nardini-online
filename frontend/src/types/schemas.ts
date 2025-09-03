// Mirror of Python schemas.py for TypeScript
export type SequenceStatus = "pending" | "cached" | "pending_external" | "complete";
export type RunStatus = "pending" | "complete";

export interface SequenceData {
  sequence_id: string;
  status: SequenceStatus;
  start_time?: number;
  end_time?: number;
  seq_uuid?: string;
  zip_path?: string;
  job_id?: string;
}

export interface RunData {
  status: RunStatus;
  fasta_filename: string;
  output_filename: string;
  sequences: Record<string, SequenceData>;
  total_sequences: number;
  cached_sequences: number;
  merged_zip_filename?: string;
  submitted_at: number;
  completed_at?: number;
}

export interface ApiResponse<T> {
  runs?: T[];
  run?: T;
  error?: string;
}