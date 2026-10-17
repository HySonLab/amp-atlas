/* Pure helical-wheel model. Sequence positions, not angular sorting, define N→C order. */
(function(scope){
'use strict';
const AA='ACDEFGHIKLMNPQRSTVWY',MAX=72;
function length(n){n=Number(n);if(!Number.isInteger(n)||n<1||n>MAX)throw Error(`Choose a whole number from 1 to ${MAX} residues.`);return n}
function validate(slots){if(!Array.isArray(slots))throw Error('Wheel positions must be an array.');length(slots.length);return slots.map(a=>{if(a===''||a===null)return '';if(typeof a!=='string'||a.length!==1||!AA.includes(a))throw Error('Use one canonical amino acid per position.');return a})}
function create(n){return Array(length(n)).fill('')}
function resize(slots,n){const s=validate(slots);return Array.from({length:length(n)},(_,i)=>s[i]||'')}
function set(slots,index,residue){const s=validate(slots);if(!Number.isInteger(index)||index<0||index>=s.length)throw Error('Invalid wheel position.');s[index]=residue;return validate(s)}
function fromSequence(text){if(typeof text!=='string')throw Error('Sequence must be text.');const s=text.replace(/\s/g,'').toUpperCase();length(s.length);if(!/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(s))throw Error('Import 1–72 canonical amino acids; FASTA headers and nonstandard symbols are not supported.');return [...s]}
function pattern(slots){return validate(slots).map(a=>a||'?').join('')}
function toSequence(slots){const s=validate(slots);if(s.some(a=>!a))throw Error('Fill every position before converting to a complete sequence.');return s.join('')}
function geometry(n){n=length(n);const outer=105+Math.floor((n-1)/18)*58,size=2*(outer+38),center=size/2;return {size,center,outer,positions:Array.from({length:n},(_,i)=>{const angle=(i*100-90)*Math.PI/180,radius=105+Math.floor(i/18)*58;return {index:i,number:i+1,ring:Math.floor(i/18),angleDegrees:(i*100)%360,x:center+radius*Math.cos(angle),y:center+radius*Math.sin(angle)}})}}
function matches(sequence,slots){const s=validate(slots);return typeof sequence==='string'&&sequence.length===s.length&&/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(sequence)&&s.every((a,i)=>!a||sequence[i]===a)}
scope.AMPWheel={AA,MAX,length,validate,create,resize,set,fromSequence,pattern,toSequence,geometry,matches};
})(typeof self!=='undefined'?self:globalThis);
