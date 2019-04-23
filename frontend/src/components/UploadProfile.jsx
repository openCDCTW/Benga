import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import blue from '@material-ui/core/colors/blue';
//
import Radio from '@material-ui/core/Radio';
import SearchBar from './SearchBar.jsx';
import ReplyIcon from '@material-ui/icons/Reply';

const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(blue[600]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
    },
    radio:{
        '&$checked': {
            color: blue[600],
    },
},

})

class Upload_profile extends React.Component {

    constructor(props) {

        super(props);

        this.state = {
            algorithm: 'single'
        };

        this.djsConfig = {
            dictDefaultMessage:"Drag cgMLST profiles here (up to 500 files)",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 500,
            init:function(){
                this.on("addedfile", function(on_load_header_data){

                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.clusteringID);
                });
                this.on("success", function(file){
                    file._removeLink.remove();
                    delete file._removeLink;
                });

            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/dendrogram/profile/'
        };

        this.dropzone = null;
    }

    handlePost() {

        let fileCheck = this.dropzone.files.length;

        if(fileCheck < 1){
            alert('Please upload at least 1 files');
            return ;
        }else if(fileCheck > 500){
            alert('Cannot upload more than 500 files');
            return ;
        }

        fetch('api/dendrogram/upload/', {method:'POST'})
            .then(function(res){
               return res.json();
            }).then(batch => window.clusteringID = batch.id);

        function submit(){

            if( window.clusteringID != undefined ){
                this.dropzone.processQueue();

                let scheme = {};
                scheme.prof_num = this.dropzone.files.length;
                scheme.linkage = this.state.algorithm;
                fetch('api/dendrogram/upload/' + window.clusteringID + '/', { 
                    method:'PATCH',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });

                this.props.history.push("/dendrogram_view");
                clearInterval(interval);
            }
        };

        let interval = setInterval(submit.bind(this),500);
    }

    remove(){
        this.dropzone.removeAllFiles();
    }

    handleChange(event){
        this.setState({ algorithm: event.target.value });
    }

    example(){

        fetch('api/dendrogram/upload/', {method:'POST'})
            .then(function(res){
               return res.json();
            }).then(batch => window.clusteringID = batch.id);

        function submit(){

            if( window.clusteringID != undefined ){
                let exampleFile = [
                    { name:"V.cholerae_01.tsv" },
                    { name:"V.cholerae_02.tsv" },
                    { name:"V.cholerae_03.tsv" },
                    { name:"V.cholerae_04.tsv" },
                    { name:"V.cholerae_05.tsv" },
                    { name:"V.cholerae_06.tsv" },
                    { name:"V.cholerae_07.tsv" },
                    { name:"V.cholerae_08.tsv" },
                    { name:"V.cholerae_09.tsv" },
                    { name:"V.cholerae_10.tsv" },
                    { name:"V.cholerae_11.tsv" },
                    { name:"V.cholerae_12.tsv" },
                ];

                var scheme = {};
                scheme.prof_num = 12;
                scheme.linkage = this.state.algorithm;
                fetch('api/dendrogram/upload/' + window.clusteringID + '/', { 
                    method:'PATCH',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });

                let i = 0;
                for(i; i < exampleFile.length; i++){
                    this.dropzone.emit("addedfile", exampleFile[i]);
                    this.dropzone.emit("success", exampleFile[i]);
                    this.dropzone.emit("complete", exampleFile[i]);
                    this.dropzone.files.push(exampleFile[i]);
                    window.fileName.push(exampleFile[i].name);
                };

                let encodeExampleData = [
                    require('./static/Example_data/V.cholerae_01.tsv'), 
                    require('./static/Example_data/V.cholerae_02.tsv'), 
                    require('./static/Example_data/V.cholerae_03.tsv'), 
                    require('./static/Example_data/V.cholerae_04.tsv'), 
                    require('./static/Example_data/V.cholerae_05.tsv'), 
                    require('./static/Example_data/V.cholerae_06.tsv'),
                    require('./static/Example_data/V.cholerae_07.tsv'), 
                    require('./static/Example_data/V.cholerae_08.tsv'), 
                    require('./static/Example_data/V.cholerae_09.tsv'), 
                    require('./static/Example_data/V.cholerae_10.tsv'), 
                    require('./static/Example_data/V.cholerae_11.tsv'), 
                    require('./static/Example_data/V.cholerae_12.tsv'), 

                ];

                let decodeExampleData = [];
                let j = 0;

                for(j; j < encodeExampleData.length; j++){
                    let tmp = encodeExampleData[j].substring(38,encodeExampleData[j].length);
                    tmp = window.atob(tmp);
                    decodeExampleData.push(tmp);
                };

                let VC01 = new File([decodeExampleData[0]],'V.cholerae_01.tsv');
                let VC02 = new File([decodeExampleData[1]],'V.cholerae_02.tsv');
                let VC03 = new File([decodeExampleData[2]],'V.cholerae_03.tsv');
                let VC04 = new File([decodeExampleData[3]],'V.cholerae_04.tsv');
                let VC05 = new File([decodeExampleData[4]],'V.cholerae_05.tsv');
                let VC06 = new File([decodeExampleData[5]],'V.cholerae_06.tsv');
                let VC07 = new File([decodeExampleData[6]],'V.cholerae_07.tsv');
                let VC08 = new File([decodeExampleData[7]],'V.cholerae_08.tsv');
                let VC09 = new File([decodeExampleData[8]],'V.cholerae_09.tsv');
                let VC10 = new File([decodeExampleData[9]],'V.cholerae_10.tsv');
                let VC11 = new File([decodeExampleData[10]],'V.cholerae_11.tsv');
                let VC12 = new File([decodeExampleData[11]],'V.cholerae_12.tsv');

                let decodeExampleFile = [ 
                    VC01, VC02, VC03, VC04, VC05, VC06, 
                    VC07, VC08, VC09, VC10, VC11, VC12
                ];

                let k = 0;
                let upload_status = [];

                for(k; k < decodeExampleFile.length; k++){
                    let form = new FormData();
                    form.append('file',decodeExampleFile[k]);
                    form.append('batch_id', window.clusteringID);

                    fetch('api/dendrogram/profile/', {
                        method:'POST',
                        body:form ,
                    });
                };

                this.props.history.push("/dendrogram_view");
                clearInterval(interval);

                }
            };

            let interval = setInterval(submit.bind(this),500);
    }

    query(){

        fetch('api/dendrogram/dendrogram/' + this.state.querybyID + '/', { method:'GET'})
        .then(function(response){
            if(response.status != 404){
                return response.json();
            }else{
                return response.status;
            }
        }).then(res => this.setState(state => ({ tmp: res })));

        function result(){

            if(this.state.tmp.png_file != undefined){
                this.setState(state => ({
                png_file: this.state.tmp.png_file, 
                pdf_file: this.state.tmp.pdf_file,
                svg_file: this.state.tmp.svg_file, 
                emf_file: this.state.tmp.emf_file, 
                newick_file: this.state.tmp.newick_file }))
            }else{
                alert("Data not found. Please input correct ID or try again later.");
            }
            clearInterval(interval);
        }
        let interval = setInterval(result.bind(this),500);
    }

    back(){
        this.setState(state => ({ png_file: undefined }));
    }

    render() {
        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,
        }
        const { classes } = this.props;

        if(this.state.png_file == undefined){
            return (
                <div>
                    <div>
                        <div style={{ float:'right', marginTop:'35px', marginRight:'25px' }}>
                            <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                                Remove all files
                                &nbsp;&nbsp;
                                <DeleteIcon />
                            </Button>
                        </div>
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <br />
                    <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                        <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                        <font>Algorithms : </font>
                        &nbsp;&nbsp;
                        <Radio
                            color='primary'
                            checked={this.state.algorithm === 'single'}
                            onChange={this.handleChange.bind(this)}
                            value="single"
                        />
                        <font>Single linkage</font>
                        <Radio
                            color='primary'
                            checked={this.state.algorithm === 'average'}
                            onChange={this.handleChange.bind(this)}
                            value="average"
                        />
                        <font>UPGMA</font>
                    </div>
                    <br />
                    <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                        <Button variant="contained" color="default" 
                         onClick={this.example.bind(this)}>
                            Example
                            &nbsp;&nbsp;
                            <CloudUploadIcon />
                        </Button>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" className ={classes.cssRoot}
                        onClick={this.handlePost.bind(this)}>
                            Submit
                            &nbsp;&nbsp;
                            <CloudUploadIcon />
                        </Button>
                    </div>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <SearchBar
                            onChange = {(value) => this.setState({ querybyID: value })}
                            onRequestSearch={this.query.bind(this)}
                            placeholder = {"Input job ID to get result"}
                            style = {{
                                width: '90%',
                                margin: '0 auto',
                            }}
                        />
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <br />
                </div>
            );
        }else{
            return (
                <div>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <img src={this.state.svg_file} />
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font>Download</font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Png</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Pdf</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">Svg</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">emf</Button>
                        </a>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">newick</Button>
                        </a>
                    </div>
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <Button variant="contained" color="default" onClick={this.back.bind(this)}>
                                <ReplyIcon />
                                &nbsp;&nbsp;
                                Back
                            </Button>
                    </div>
                </div>
            );
        }
    }
}

export default withStyles(styles)(Upload_profile);