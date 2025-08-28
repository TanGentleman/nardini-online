# Nardini Online - Technical Guide

This guide provides comprehensive technical documentation for deploying, configuring, and using the Nardini Online API.

## Table of Contents

- [Architecture Overview](#architecture-overview)
- [Installation & Deployment](#installation--deployment)
- [API Reference](#api-reference)
- [Configuration](#configuration)
- [Development Setup](#development-setup)
- [Troubleshooting](#troubleshooting)

## Architecture Overview

Nardini Online is built on **Modal**, a serverless computing platform, providing scalable protein sequence analysis using the NARDINI tool for Intrinsically Disordered Regions (IDRs).

### Core Components

- **FastAPI Web Service**: Handles file uploads, job management, and result downloads
- **Distributed Processing**: Sequences processed in parallel using Modal's serverless functions
- **Persistent Storage**: Modal Volume for caching results and storing metadata
- **Batch Processing**: Optimized batching (16 sequences per worker) for efficient resource usage

### Processing Flow

1. **Upload**: FASTA file uploaded via REST API
2. **Parsing**: Sequences extracted and validated
3. **Caching**: Check existing results to avoid reprocessing
4. **Processing**: Novel sequences submitted for parallel analysis
5. **Storage**: Results cached as individual zip files
6. **Merging**: Final results merged into single downloadable archive

## Installation & Deployment

### Prerequisites

- Python 3.9+
- Modal account and CLI setup
- Required dependencies (see `requirements.txt`)

### Environment Variables

Configure these environment variables for customization:

```bash
APP_NAME=nardini_online           # Modal app name
VOLUME_NAME=run_fasta_volume      # Modal volume name
TIMEOUT_SECONDS=21600             # Processing timeout (6 hours default)
MAX_UPLOAD_MB=10                  # Maximum file upload size
```

### Deployment Steps

1. **Install Modal CLI**:
   ```bash
   pip install modal
   modal setup
   ```

2. **Clone Repository**:
   ```bash
   git clone <repository-url>
   cd nardini-online
   ```

3. **Deploy to Modal**:
   ```bash
   modal deploy src/app/backend.py
   ```

4. **Get API URL**:
   After deployment, Modal provides the FastAPI endpoint URL.

## API Reference

Base URL: `https://your-app-name--fastapi-app.modal.run`

### Endpoints

#### Health Check
```http
GET /health
```

**Response**:
```json
{
  "status": "healthy"
}
```

#### Upload FASTA File
```http
POST /upload_fasta
Content-Type: multipart/form-data
```

**Parameters**:
- `file` (required): FASTA file (.fasta, .fa, .fas)
- `output_filename` (optional): Custom output filename

**Response**:
```json
{
  "run_id": "uuid-string",
  "status": "submitted|ready",
  "message": "Success message",
  "job_ids": ["job1", "job2", "..."]
}
```

**Status Values**:
- `submitted`: Novel sequences are being processed
- `ready`: All sequences were cached (no processing needed)

#### Check Job Status
```http
GET /status/{run_id}
```

**Response**:
```json
{
  "run_id": "uuid-string",
  "status": "pending|complete",
  "pending_sequences": ["seq_id1", "seq_id2"]
}
```

#### Download Results
```http
GET /download/{run_id}?output_filename=custom_name.zip
```

**Parameters**:
- `output_filename` (optional): Override the default output filename

**Response**: ZIP file containing:
- Individual TSV files with statistical analysis
- PNG visualization files
- Merged summary data

#### Retry Failed Sequences
```http
GET /retry/{run_id}
```

**Response**:
```json
{
  "run_id": "uuid-string",
  "status": "retry_submitted"
}
```

### Error Responses

All endpoints return error responses in this format:

```json
{
  "detail": "Error message description"
}
```

Common HTTP status codes:
- `400`: Bad request (invalid file, missing parameters)
- `404`: Run not found
- `413`: File too large
- `500`: Internal server error

## Configuration

### File Size Limits

- Maximum upload: 10MB (configurable via `MAX_UPLOAD_MB`)
- Supported formats: `.fasta`, `.fa`, `.fas`

### Processing Limits

- Timeout: 6 hours per job (configurable via `TIMEOUT_SECONDS`)
- Batch size: 16 sequences per worker (optimized for performance)
- Concurrent workers: Unlimited (Modal's auto-scaling)

### Caching Behavior

- Results cached by sequence string hash
- Duplicate sequences across runs reuse cached results
- Cache persists across deployments (stored in Modal Volume)

## Development Setup

### Local Development

1. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run Local Tests**:
   ```bash
   modal run src/app/backend.py::test_function
   ```

3. **Debug Mode**:
   Set `REQUIRE_AUTH = True` in `backend.py` for authentication during development.

### Project Structure

```
src/
├── app/
│   ├── backend.py          # Main FastAPI application
│   └── common.py           # Shared Modal configuration
├── shared_utils/
│   ├── schemas.py          # TypedDict definitions
│   ├── file_utils.py       # File system operations
│   └── utils.py            # Utility functions
├── backend_utils/
│   └── utils.py            # Processing logic
└── fastapi_utils/
    └── utils.py            # API helper functions
```

### Key Functions

- `process_single_sequence()`: Processes individual sequences
- `process_16_sequences()`: Batch processing for efficiency  
- `retry_pending_sequences()`: Handles failed job recovery
- `migrate_run_metadata()`: Updates metadata format

## Troubleshooting

### Common Issues

#### 1. "Run not found" Error
- Verify the `run_id` is correct
- Check if the run metadata exists in the Modal Volume

#### 2. Processing Stuck in "Pending"
- Use `/retry/{run_id}` endpoint to restart failed jobs
- Check Modal dashboard for worker errors
- Verify sequence format is valid

#### 3. Download Fails
- Ensure all sequences are complete before downloading
- Check `/status/{run_id}` first to verify completion
- Verify sufficient storage space in Modal Volume

#### 4. Upload Rejected
- Check file size (max 10MB default)
- Verify file extension (.fasta, .fa, .fas)
- Ensure file contains valid FASTA sequences

### Monitoring

#### Modal Dashboard
- Monitor function invocations and errors
- View logs for debugging
- Check volume usage and costs

#### Logging
- All operations logged with structured messages
- Error details captured for troubleshooting
- Performance metrics tracked

### Performance Optimization

#### Batch Size Tuning
- Default 16 sequences per worker optimized for most use cases
- Adjust `BATCH_SIZE` in code for specific workloads
- Consider sequence length distribution

#### Caching Strategy
- Sequences cached by exact string match
- Consider sequence normalization for better hit rates
- Monitor cache hit/miss ratios

## Security Considerations

- Files validated for FASTA format before processing
- No sensitive data logged or exposed
- Rate limiting available through Modal's infrastructure
- Optional authentication via `REQUIRE_AUTH` flag

## Example Usage

### Python Client Example

```python
import requests
import time

# API base URL (replace with your actual deployment URL)
BASE_URL = "https://your-app-name--fastapi-app.modal.run"

def process_fasta_file(fasta_file_path, output_filename=None):
    """Complete workflow example for processing a FASTA file."""
    
    # 1. Health check
    response = requests.get(f"{BASE_URL}/health")
    print(f"API Status: {response.json()['status']}")
    
    # 2. Upload FASTA file
    with open(fasta_file_path, 'rb') as f:
        files = {'file': f}
        params = {}
        if output_filename:
            params['output_filename'] = output_filename
            
        response = requests.post(f"{BASE_URL}/upload_fasta", files=files, params=params)
    
    if response.status_code != 200:
        print(f"Upload failed: {response.json()}")
        return None
        
    upload_result = response.json()
    run_id = upload_result['run_id']
    print(f"Upload successful. Run ID: {run_id}")
    print(f"Status: {upload_result['status']}")
    
    # 3. Poll for completion
    while True:
        response = requests.get(f"{BASE_URL}/status/{run_id}")
        status_result = response.json()
        
        if status_result['status'] == 'complete':
            print("Processing complete!")
            break
        elif status_result['status'] == 'pending':
            pending_count = len(status_result['pending_sequences'])
            print(f"Still processing {pending_count} sequences...")
            time.sleep(10)  # Wait 10 seconds before checking again
        
    # 4. Download results
    response = requests.get(f"{BASE_URL}/download/{run_id}")
    if response.status_code == 200:
        output_file = output_filename or f"results_{run_id}.zip"
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print(f"Results downloaded to: {output_file}")
        return output_file
    else:
        print(f"Download failed: {response.status_code}")
        return None

# Usage
if __name__ == "__main__":
    result_file = process_fasta_file("my_sequences.fasta", "my_results.zip")
    if result_file:
        print(f"Processing complete! Results saved to {result_file}")
```

### cURL Examples

#### Upload FASTA File
```bash
curl -X POST "https://your-app-name--fastapi-app.modal.run/upload_fasta" \
  -H "Content-Type: multipart/form-data" \
  -F "file=@sequences.fasta" \
  -F "output_filename=my_results.zip"
```

#### Check Status
```bash
curl "https://your-app-name--fastapi-app.modal.run/status/your-run-id-here"
```

#### Download Results
```bash
curl -O "https://your-app-name--fastapi-app.modal.run/download/your-run-id-here?output_filename=results.zip"
```

### JavaScript/Node.js Example

```javascript
const axios = require('axios');
const fs = require('fs');
const FormData = require('form-data');

const BASE_URL = 'https://your-app-name--fastapi-app.modal.run';

async function processFasta(fastaFilePath, outputFilename) {
    try {
        // 1. Upload file
        const form = new FormData();
        form.append('file', fs.createReadStream(fastaFilePath));
        if (outputFilename) {
            form.append('output_filename', outputFilename);
        }
        
        const uploadResponse = await axios.post(`${BASE_URL}/upload_fasta`, form, {
            headers: form.getHeaders()
        });
        
        const runId = uploadResponse.data.run_id;
        console.log(`Upload successful. Run ID: ${runId}`);
        
        // 2. Poll for completion
        while (true) {
            const statusResponse = await axios.get(`${BASE_URL}/status/${runId}`);
            
            if (statusResponse.data.status === 'complete') {
                console.log('Processing complete!');
                break;
            }
            
            console.log(`Still processing ${statusResponse.data.pending_sequences.length} sequences...`);
            await new Promise(resolve => setTimeout(resolve, 10000)); // Wait 10 seconds
        }
        
        // 3. Download results
        const downloadResponse = await axios.get(`${BASE_URL}/download/${runId}`, {
            responseType: 'stream'
        });
        
        const outputPath = outputFilename || `results_${runId}.zip`;
        const writer = fs.createWriteStream(outputPath);
        downloadResponse.data.pipe(writer);
        
        return new Promise((resolve, reject) => {
            writer.on('finish', () => {
                console.log(`Results downloaded to: ${outputPath}`);
                resolve(outputPath);
            });
            writer.on('error', reject);
        });
        
    } catch (error) {
        console.error('Error:', error.response?.data || error.message);
        throw error;
    }
}

// Usage (best if sequences already cached)
processFasta('sequences.fasta', 'my_results.zip')
    .then(outputPath => console.log(`Processing complete! Results: ${outputPath}`))
    .catch(error => console.error('Failed:', error));
```

## Support

For issues or questions:

1. Check the [troubleshooting section](#troubleshooting)
2. Review Modal platform documentation
3. Contact the development team

---

**Built with ❤️ by Tanuj Vasudeva and Ethan Caine, 2025**
